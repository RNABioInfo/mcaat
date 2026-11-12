#ifndef MCAAT_REACHABLE_MEMBERSHIP_H
#define MCAAT_REACHABLE_MEMBERSHIP_H
#include "protospacer_detection/extension_engine.h"
namespace mcaat::proto::audit {
// Sound relaxation of membership: reject only nodes proved disconnected from
// every start. A bounded reverse search that runs out of work returns UNKNOWN
// (accepted), never a negative. This is not shortest-path sequence recovery.
template<class Graph> struct ReachableMembership {
    Graph& graph;
    std::set<Node> roots,disconnected,accepted;
    bool backward;
    size_t work_limit=2000,queries=0,negative=0,unknown=0;
    ReachableMembership(Graph& g,const Path& r,bool back):graph(g),roots(r.begin(),r.end()),backward(back){}
    size_t k()const{return graph.k();}
    Node find(const std::string& word) {
        auto node=graph.find(word);if(node==missing || disconnected.count(node))return missing;
        if(roots.count(node) || accepted.count(node))return node;
        ++queries;
        std::set<Node> seen{node};Path queue{node};
        for(size_t i=0;i<queue.size();++i) {
            for(auto arc:graph.neighbors(queue[i],!backward)) {
                if(roots.count(arc.to) || accepted.count(arc.to)){accepted.insert(node);return node;}
                if(disconnected.count(arc.to) || !seen.insert(arc.to).second)continue;
                queue.push_back(arc.to);
                if(queue.size()>=work_limit){++unknown;accepted.insert(node);return node;}
            }
        }
        ++negative;disconnected.insert(seen.begin(),seen.end());return missing;
    }
};
}
#endif
