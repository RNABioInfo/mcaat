#ifndef MCAAT_REFERENCE_GRAPH_ORACLE_H
#define MCAAT_REFERENCE_GRAPH_ORACLE_H
// Diagnostic only. References must never be passed to the production ranker.
#include "protospacer_detection/extension_engine.h"
#include <queue>
namespace mcaat::proto::audit {
struct Key {Node node;size_t pos;bool operator==(const Key& b)const{return node==b.node&&pos==b.pos;}};
struct KeyHash {size_t operator()(const Key& k)const{return std::hash<Node>{}(k.node)^(k.pos*0x9e3779b97f4a7c15ULL);}};
struct State {Node node;size_t pos,cost,parent;char op;};
struct Witness {Path nodes;std::string sequence,operations;size_t ref_bp=0,edits=0,matches=0,columns=0;};
struct Search {
    std::map<size_t,Witness> found;
    size_t popped=0,discovered=0;
    size_t lowest_unexplored_cost=std::numeric_limits<size_t>::max();
    bool limited=false,exhausted=false;
};
template<class Graph> Search align(Graph& g,Node start,const std::string& reference,bool backward,
                                  const std::set<size_t>& targets,size_t max_states) {
    Search out;if(reference.empty()||start==missing)return out;
    // Dijkstra on (graph vertex, reference position). Match=0, substitution/deletion=10, insertion=9.
    // At reference length L, cost<=L is exactly >=90% aligned identity:
    // matches/(L+insertions)>=0.9. Revisits are permitted; (position,cost) bounds
    // the finite search. No multiplicity gate, source filter or live beam.
    size_t budget=reference.size();
    std::vector<State> states{{start,0,0,size_t(-1),'S'}};
    std::unordered_map<Key,size_t,KeyHash> best;best.emplace(Key{start,0},0);
    using Item=std::tuple<size_t,size_t,size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    queue.emplace(0,reference.size(),0);
    std::unordered_map<Node,std::vector<std::pair<Arc,char>>> adjacency;
    auto push=[&](Node v,size_t pos,size_t cost,size_t parent,char op) {
        if(cost>budget || pos>reference.size())return;
        Key key{v,pos};auto it=best.find(key);
        if(it!=best.end()&&states[it->second].cost<=cost)return;
        if(states.size()>=max_states){out.limited=true;out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,cost);return;}
        size_t index=states.size();states.push_back({v,pos,cost,parent,op});best[key]=index;
        queue.emplace(cost,reference.size()-pos,index);
    };
    while(!queue.empty()) {
        auto [cost,remaining,index]=queue.top();queue.pop();auto state=states[index];
        if(best.at(Key{state.node,state.pos})!=index)continue;
        ++out.popped;
        if(targets.count(state.pos) && cost<=state.pos && !out.found.count(state.pos)) {
            std::vector<size_t> chain;for(size_t i=index;states[i].parent!=size_t(-1);i=states[i].parent)chain.push_back(i);
            std::reverse(chain.begin(),chain.end());Witness w;w.nodes={start};w.ref_bp=state.pos;w.edits=0;
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
        }
        bool complete=true;for(auto t:targets)if(!out.found.count(t))complete=false;
        if(complete)break;
        if(state.pos==reference.size())continue;
        // A state already over every remaining target's edit budget cannot
        // contribute a qualifying witness. This is a sound cost bound.
        bool useful=false;for(auto t:targets)if(t>=state.pos&&!out.found.count(t)&&cost<=t){useful=true;break;}
        if(!useful)continue;
        auto it=adjacency.find(state.node);
        if(it==adjacency.end()) {
            std::vector<std::pair<Arc,char>> choices;
            for(auto a:g.neighbors(state.node,backward))choices.push_back({a,g.extension_base(a,backward)});
            it=adjacency.emplace(state.node,std::move(choices)).first;
        }
        push(state.node,state.pos+1,cost+10,index,'D');
        for(const auto& choice:it->second) {
            const auto& arc=choice.first;bool same=choice.second==reference[state.pos];
            push(arc.to,state.pos+1,cost+(same?0:10),index,same?'M':'X');
            push(arc.to,state.pos,cost+9,index,'I');
        }
    }
    if(!queue.empty())out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,std::get<0>(queue.top()));
    out.discovered=states.size();out.exhausted=queue.empty()&&!out.limited;return out;
}
} // namespace mcaat::proto::audit
#endif
