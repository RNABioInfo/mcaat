#ifndef MCAAT_REFERENCE_CYCLE_CONTEXT_H
#define MCAAT_REFERENCE_CYCLE_CONTEXT_H
#include "protospacer_detection/extension_engine.h"
namespace mcaat::proto::audit {
// Exact state for the source-edge closed-walk exclusion. The sorted set stores
// previous vertices since the last non-source edge, with the current vertex
// supplied by the alignment state. Other repeats are unrestricted here.
class CycleContext {
    SourceCycleGuard guard_;
    std::vector<Path> previous_{{}};
    std::map<Path,size_t> ids_{{{},0}};
public:
    static constexpr size_t reject=missing-1;
    CycleContext(size_t k,const std::vector<FlankContext>& models):guard_(k,&models){}
    template<class Graph> size_t advance(Graph& graph,size_t state,Node from,const Arc& arc,bool back) {
        if(!guard_.contains(graph,from,arc,back))return 0;
        const auto& before=previous_.at(state);
        if(arc.to==from || std::binary_search(before.begin(),before.end(),arc.to))return reject;
        auto next=before;next.insert(std::lower_bound(next.begin(),next.end(),from),from);
        auto found=ids_.find(next);if(found!=ids_.end())return found->second;
        size_t id=previous_.size();ids_.emplace(next,id);previous_.push_back(std::move(next));return id;
    }
};
}
#endif
