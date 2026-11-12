#ifndef MCAAT_REFERENCE_UNITIG_ORACLE_H
#define MCAAT_REFERENCE_UNITIG_ORACLE_H
#include "reference_graph_oracle.h"
#include <memory>
namespace mcaat::proto::audit {
// Exact graph/reference alignment with linear stretches handled by string DP.
// Every reference endpoint remains eligible inside a stretch. Junction states
// keep reference position and both path-dependent source-context states.
template<class Graph> Search align_unitigs(Graph& graph,Node start,const std::string& reference,bool back,
        const std::set<size_t>& targets,size_t max_states,
        const std::function<size_t(size_t,size_t)>& suffix_bound={},const Path& extra_starts={},
        bool first_only=false,bool prefer_progress=false,PrefixContext* context=nullptr,CycleContext* cycles=nullptr) {
    Search out;if(reference.empty() || (start==missing&&extra_starts.empty()))return out;
    struct Span {Path nodes;std::string bases;std::vector<Arc> arcs;};
    struct Vertex {Node node;size_t pos,cost,parent,context,cycle;std::shared_ptr<Span> span;size_t used=0;};
    const size_t inf=std::numeric_limits<size_t>::max()/4,budget=reference.size();
    std::vector<int64_t> admissible(reference.size()+1,-1);
    for(size_t pos=0;pos<=reference.size();++pos)for(auto target:targets)if(target>=pos)
        admissible[pos]=std::max(admissible[pos],int64_t(target)-int64_t(suffix_bound?suffix_bound(pos,target):0));
    std::vector<Vertex> states;std::vector<bool> settled;std::unordered_map<Key,size_t,KeyHash> best;
    using Item=std::tuple<int64_t,size_t,size_t>;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;
    auto push=[&](Vertex state) {
        if(int64_t(state.cost)>admissible[state.pos])return;
        Key key{state.node,state.pos,state.context,state.cycle};auto it=best.find(key);
        if(it!=best.end() && states[it->second].cost<=state.cost)return;
        // In Dijkstra order, no child can refer to an unpopped state. Updating
        // it in place is safe and avoids charging the work/memory cap for each
        // obsolete proposal to the same junction/reference-position key.
        if(!prefer_progress && it!=best.end()) {
            auto index=it->second;if(settled[index])throw std::runtime_error("negative-cost unitig relaxation");
            states[index]=std::move(state);const auto& saved=states[index];
            queue.emplace(saved.cost,reference.size()-saved.pos,index);return;
        }
        if(states.size()>=max_states){out.limited=true;out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,state.cost);return;}
        size_t index=states.size();best[key]=index;states.push_back(std::move(state));settled.push_back(false);const auto& saved=states.back();
        queue.emplace(int64_t(saved.cost)-(prefer_progress?5*int64_t(saved.pos):0),reference.size()-saved.pos,index);
    };
    std::set<Node> roots(extra_starts.begin(),extra_starts.end());if(start!=missing)roots.insert(start);
    for(auto root:roots)push({root,0,0,size_t(-1),context?context->initial():0,0,{},0});
    std::map<std::pair<Node,Node>,std::shared_ptr<Span>> spans;
    auto get_span=[&](Node from,Arc first) {
        auto key=std::make_pair(from,first.edge);auto cached=spans.find(key);if(cached!=spans.end())return cached->second;
        auto span=std::make_shared<Span>();span->nodes.push_back(from);std::set<Node> visited{from};auto arc=first;
        while(true) {
            span->arcs.push_back(arc);span->nodes.push_back(arc.to);span->bases+=graph.extension_base(arc,back);
            if(!visited.insert(arc.to).second || span->bases.size()>=reference.size()+budget/9)break;
            auto next=graph.neighbors(arc.to,back);
            if(next.size()!=1 || graph.neighbors(arc.to,!back).size()!=1)break;
            arc=next.front();
        }
        spans.emplace(key,span);return span;
    };
    auto reconstruct=[&](const Vertex& terminal) {
        std::vector<Vertex> chain{terminal};auto parent=terminal.parent;
        while(parent!=size_t(-1)){chain.push_back(states[parent]);parent=states[parent].parent;}
        std::reverse(chain.begin(),chain.end());Witness witness;witness.nodes={chain.front().node};witness.ref_bp=terminal.pos;
        size_t total_cost=0;
        for(size_t segment=1;segment<chain.size();++segment) {
            const auto& current=chain[segment];const auto& prior=chain[segment-1];
            std::string seq=current.span?current.span->bases.substr(0,current.used):"";
            auto ref=reference.substr(prior.pos,current.pos-prior.pos);size_t width=ref.size()+1;
            std::vector<size_t> matrix((seq.size()+1)*width);
            auto at=[&](size_t i,size_t j)->size_t&{return matrix[i*width+j];};
            for(size_t j=0;j<=ref.size();++j)at(0,j)=10*j;
            for(size_t i=1;i<=seq.size();++i){at(i,0)=9*i;for(size_t j=1;j<=ref.size();++j)
                at(i,j)=std::min({at(i-1,j)+9,at(i,j-1)+10,at(i-1,j-1)+(seq[i-1]==ref[j-1]?0:10)});}
            total_cost+=at(seq.size(),ref.size());
            size_t i=seq.size(),j=ref.size();std::string ops;
            while(i||j) {
                if(i&&j&&at(i,j)==at(i-1,j-1)+(seq[i-1]==ref[j-1]?0:10)){ops+=seq[i-1]==ref[j-1]?'M':'X';--i;--j;}
                else if(i&&at(i,j)==at(i-1,j)+9){ops+='I';--i;}
                else{ops+='D';--j;}
            }
            std::reverse(ops.begin(),ops.end());witness.operations+=ops;witness.sequence+=seq;
            if(current.span)witness.nodes.insert(witness.nodes.end(),current.span->nodes.begin()+1,current.span->nodes.begin()+1+current.used);
        }
        for(char op:witness.operations){++witness.columns;if(op=='M')++witness.matches;else ++witness.edits;}
        if(total_cost!=terminal.cost || witness.matches*100<90*witness.columns)throw std::runtime_error("unitig oracle reconstructed an invalid alignment");
        auto path=witness.nodes;if(back)std::reverse(path.begin(),path.end());auto spelled=graph.spell(path);
        auto emitted=back?spelled.substr(0,spelled.size()-graph.k()):spelled.substr(graph.k());if(back)std::reverse(emitted.begin(),emitted.end());
        if(emitted!=witness.sequence)throw std::runtime_error("unitig oracle path does not spell its witness");
        return witness;
    };
    auto accept=[&](const Vertex& state) {
        if(state.context || !targets.count(state.pos) || state.cost>state.pos || out.found.count(state.pos))return false;
        out.found.emplace(state.pos,reconstruct(state));
        return first_only || out.found.size()==targets.size();
    };
    bool complete=false;
    while(!queue.empty()&&!complete) {
        auto [priority,remaining,index]=queue.top();queue.pop();auto state=states[index];
        if(best.at(Key{state.node,state.pos,state.context,state.cycle})!=index ||
            (!prefer_progress && (settled[index] || priority!=int64_t(state.cost))))continue;
        settled[index]=true;++out.popped;
        if(accept(state)){complete=true;break;}
        // A walk can end at this vertex while deleting the remaining reference
        // prefix. No graph transition is necessary for such an endpoint.
        for(auto target:targets)if(target>state.pos && state.cost+10*(target-state.pos)<=target) {
            Vertex end{state.node,target,state.cost+10*(target-state.pos),index,state.context,state.cycle,{},0};
            if(accept(end)){complete=true;break;}
        }
        if(complete || state.pos==reference.size())continue;
        for(auto first:graph.neighbors(state.node,back)) {
            auto span=get_span(state.node,first);size_t length=reference.size()-state.pos;
            std::vector<size_t> row(length+1,inf),next(length+1,inf);
            for(size_t j=0;j<=length;++j)if(int64_t(state.cost+10*j)<=admissible[state.pos+j])row[j]=10*j;
            auto prefix=state.context,cycle=state.cycle;bool reached_end=false;
            for(size_t i=1;i<=span->bases.size();++i) {
                if(context){prefix=context->advance(prefix,span->bases[i-1]);if(prefix==PrefixContext::reject)break;}
                if(cycles){cycle=cycles->advance(graph,cycle,span->nodes[i-1],span->arcs[i-1],back);if(cycle==CycleContext::reject)break;}
                next[0]=row[0]+9;
                if(next[0]>=inf || int64_t(state.cost+next[0])>admissible[state.pos])next[0]=inf;
                bool alive=next[0]<inf;
                for(size_t j=1;j<=length;++j) {
                    next[j]=std::min({row[j]+9,next[j-1]+10,row[j-1]+(span->bases[i-1]==reference[state.pos+j-1]?0:10)});
                    if(next[j]>=inf || int64_t(state.cost+next[j])>admissible[state.pos+j])next[j]=inf;
                    alive|=next[j]<inf;
                    if(next[j]<inf && accept({span->nodes[i],state.pos+j,state.cost+next[j],index,prefix,cycle,span,i})){complete=true;break;}
                }
                row.swap(next);if(complete||!alive)break;
                reached_end=i==span->bases.size();
            }
            if(complete)break;
            if(reached_end)for(size_t j=0;j<=length;++j)if(row[j]<inf)
                push({span->nodes.back(),state.pos+j,state.cost+row[j],index,prefix,cycle,span,span->bases.size()});
        }
    }
    if(!queue.empty())out.lowest_unexplored_cost=std::min(out.lowest_unexplored_cost,prefer_progress?size_t(0):size_t(std::get<0>(queue.top())));
    out.discovered=states.size();out.exhausted=queue.empty()&&!out.limited;return out;
}
}
#endif
