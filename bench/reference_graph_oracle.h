#ifndef MCAAT_REFERENCE_GRAPH_ORACLE_H
#define MCAAT_REFERENCE_GRAPH_ORACLE_H
// Diagnostic only. References must never be passed to the production ranker.
#include "protospacer_detection/extension_engine.h"
#include "reference_prefix_context.h"
#include "reference_cycle_context.h"
#include <queue>
namespace mcaat::proto::audit {
struct Key {Node node;size_t pos,context=0,cycle=0;bool operator==(const Key& b)const{return node==b.node&&pos==b.pos&&context==b.context&&cycle==b.cycle;}};
struct KeyHash {size_t operator()(const Key& k)const{return std::hash<Node>{}(k.node)^(k.pos*0x9e3779b97f4a7c15ULL)^(k.context*0xc2b2ae3d27d4eb4fULL)^(k.cycle*0x165667b19e3779f9ULL);}};
struct State {Node node;size_t pos,cost,parent;char op;size_t context=0,cycle=0;};
struct Witness {Path nodes;std::string sequence,operations;size_t ref_bp=0,edits=0,matches=0,columns=0;};
struct Search {
    std::map<size_t,Witness> found;
    size_t popped=0,discovered=0;
    size_t lowest_unexplored_cost=std::numeric_limits<size_t>::max();
    bool limited=false,exhausted=false;
};
template<class Graph> Search align(Graph& g,Node start,const std::string& reference,bool backward,
                                  const std::set<size_t>& targets,size_t max_states,
                                  const std::function<size_t(size_t,size_t)>& suffix_bound={},
                                  const Path& extra_starts={},bool first_only=false,bool prefer_progress=false,
                                  PrefixContext* context=nullptr,CycleContext* cycles=nullptr) {
    Search out;if(reference.empty()||(start==missing&&extra_starts.empty()))return out;
    // Dijkstra on (graph vertex, reference position). Match=0, substitution/deletion=10, insertion=9.
    // At reference length L, cost<=L is exactly >=90% aligned identity:
    // matches/(L+insertions)>=0.9. Revisits are permitted; (position,cost) bounds
    // the finite search. No multiplicity gate, source filter or live beam.
    size_t budget=reference.size();
    std::vector<int64_t> admissible(reference.size()+1,-1);
    if(first_only)for(size_t pos=0;pos<=reference.size();++pos)for(auto target:targets)if(target>=pos) {
        int64_t bound=suffix_bound?suffix_bound(pos,target):0;
        admissible[pos]=std::max(admissible[pos],int64_t(target)-bound);
    }
    std::vector<State> states;
    std::unordered_map<Key,size_t,KeyHash> best;
    using Item=std::tuple<int64_t,size_t,size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    std::set<Node> roots(extra_starts.begin(),extra_starts.end());if(start!=missing)roots.insert(start);
    for(auto root:roots) {
        if(states.size()>=max_states){out.limited=true;out.lowest_unexplored_cost=0;break;}
        size_t initial=context?context->initial():0;
        size_t id=states.size();states.push_back({root,0,0,size_t(-1),'S',initial});
        best.emplace(Key{root,0,initial},id);queue.emplace(0,reference.size(),id);
    }
    std::unordered_map<Node,std::vector<std::pair<Arc,char>>> adjacency;
    auto push=[&](Node v,size_t pos,size_t cost,size_t parent,char op,size_t prefix,size_t cycle) {
        if(cost>budget || pos>reference.size())return;
        if(first_only) {
            if(int64_t(cost)>admissible[pos])return;
        } else if(suffix_bound) {
            bool useful=false;
            for(auto target:targets)if(target>=pos && !out.found.count(target) &&
                    cost+suffix_bound(pos,target)<=target){useful=true;break;}
            if(!useful)return;
        }
        Key key{v,pos,prefix,cycle};auto it=best.find(key);
        if(it!=best.end()&&states[it->second].cost<=cost)return;
        if(states.size()>=max_states){out.limited=true;out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,cost);return;}
        size_t index=states.size();states.push_back({v,pos,cost,parent,op,prefix,cycle});best[key]=index;
        queue.emplace(int64_t(cost)-(prefer_progress?5*int64_t(pos):0),reference.size()-pos,index);
    };
    while(!queue.empty()) {
        auto [priority,remaining,index]=queue.top();queue.pop();auto state=states[index];auto cost=state.cost;
        if(best.at(Key{state.node,state.pos,state.context,state.cycle})!=index)continue;
        ++out.popped;
        if(!state.context && targets.count(state.pos) && cost<=state.pos && !out.found.count(state.pos)) {
            std::vector<size_t> chain;for(size_t i=index;states[i].parent!=size_t(-1);i=states[i].parent)chain.push_back(i);
            std::reverse(chain.begin(),chain.end());Witness w;
            auto root=index;while(states[root].parent!=size_t(-1))root=states[root].parent;
            w.nodes={states[root].node};w.ref_bp=state.pos;w.edits=0;
            for(auto i:chain) {
                const auto& v=states[i];w.operations+=v.op;++w.columns;
                if(v.op=='M')++w.matches;else ++w.edits;
                if(v.op!='D') {
                    w.nodes.push_back(v.node);auto label=g.label(v.node);w.sequence+=backward?label.front():label.back();
                }
            }
            if(w.matches*100<90*w.columns)throw std::runtime_error("oracle witness below identity threshold");
            auto path=w.nodes;if(backward)std::reverse(path.begin(),path.end());
            auto spelled=g.spell(path);auto observed=backward?spelled.substr(0,spelled.size()-g.k()):spelled.substr(g.k());
            if(backward)std::reverse(observed.begin(),observed.end());
            if(observed!=w.sequence)throw std::runtime_error("oracle witness not spelled by graph");
            out.found.emplace(state.pos,std::move(w));
            if(first_only)break;
        }
        bool complete=!first_only;
        if(!first_only)for(auto t:targets)if(!out.found.count(t))complete=false;
        if(complete)break;
        if(state.pos==reference.size())continue;
        // A state already over every remaining target's edit budget cannot
        // contribute a qualifying witness. This is a sound cost bound.
        bool useful=first_only && int64_t(cost)<=admissible[state.pos];
        if(!first_only)for(auto t:targets)if(t>=state.pos&&!out.found.count(t)&&cost<=t){useful=true;break;}
        if(!useful)continue;
        auto it=adjacency.find(state.node);
        if(it==adjacency.end()) {
            std::vector<std::pair<Arc,char>> choices;
            for(auto a:g.neighbors(state.node,backward))choices.push_back({a,g.extension_base(a,backward)});
            it=adjacency.emplace(state.node,std::move(choices)).first;
        }
        push(state.node,state.pos+1,cost+10,index,'D',state.context,state.cycle);
        for(const auto& choice:it->second) {
            const auto& arc=choice.first;bool same=choice.second==reference[state.pos];
            auto prefix=context?context->advance(state.context,choice.second):0;
            if(prefix==PrefixContext::reject)continue;
            auto cycle=cycles?cycles->advance(g,state.cycle,state.node,arc,backward):0;
            if(cycle==CycleContext::reject)continue;
            push(arc.to,state.pos+1,cost+(same?0:10),index,same?'M':'X',prefix,cycle);
            push(arc.to,state.pos,cost+9,index,'I',prefix,cycle);
        }
    }
    if(!queue.empty())out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,prefer_progress?size_t(0):size_t(std::get<0>(queue.top())));
    out.discovered=states.size();out.exhausted=queue.empty()&&!out.limited;return out;
}
} // namespace mcaat::proto::audit
#endif
