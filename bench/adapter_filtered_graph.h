#ifndef MCAAT_ADAPTER_FILTERED_GRAPH_H
#define MCAAT_ADAPTER_FILTERED_GRAPH_H
#include "protospacer_detection/extension_engine.h"
namespace mcaat::proto::audit {
// Read-only induced graph. All default adapter spans fit inside a k-mer at
// k=23, so an adapter-free emitted sequence cannot visit a rejected node.
// For a short seed, whitelist its initial contexts until the whole current
// k-mer lies in the emitted seed+flank. This is an optimistic relaxation: a
// padded context outside the emitted interval must not cause a false negative.
template<class Graph> class AdapterFilteredGraph {
    Graph& graph_;
    bool enabled_;
    AdapterGuard guard_{Options{}.adapters};
    std::unordered_map<Node,bool>& blocked_;
    std::set<Node> initial_;
public:
    AdapterFilteredGraph(Graph& graph,bool enabled,const Path& roots,size_t seed_length,bool backward,
                         std::unordered_map<Node,bool>& cache):graph_(graph),enabled_(enabled),blocked_(cache) {
        if(!enabled_)return;
        if(k()<21)throw std::runtime_error("adapter node audit requires k >= 21");
        Path layer=roots;initial_.insert(roots.begin(),roots.end());
        for(size_t depth=0;seed_length+depth<k();++depth) {
            Path next;for(auto node:layer)for(auto arc:graph_.neighbors(node,backward))
                if(initial_.insert(arc.to).second)next.push_back(arc.to);
            layer=std::move(next);
        }
    }
    size_t k()const{return graph_.k();}
    bool allowed(Node node) {
        if(node==missing)return false;
        if(!enabled_ || initial_.count(node))return true;
        auto it=blocked_.find(node);
        if(it==blocked_.end())it=blocked_.emplace(node,guard_.contains(graph_.label(node))).first;
        return !it->second;
    }
    Node find(const std::string& word){auto node=graph_.find(word);return allowed(node)?node:missing;}
    std::vector<Arc> neighbors(Node node,bool back){auto arcs=graph_.neighbors(node,back);arcs.erase(std::remove_if(arcs.begin(),arcs.end(),[&](const auto& arc){return !allowed(arc.to);}),arcs.end());return arcs;}
    std::string label(Node node)const{return graph_.label(node);}
    char extension_base(const Arc& arc,bool back)const{return graph_.extension_base(arc,back);}
    std::string spell(const Path& path)const{return graph_.spell(path);}
    double count(Node a,Node b)const{return graph_.count(a,b);}
};
}
#endif
