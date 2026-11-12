#ifndef MCAAT_REFERENCE_ROOT_GRAPH_H
#define MCAAT_REFERENCE_ROOT_GRAPH_H
#include "protospacer_detection/extension_engine.h"
namespace mcaat::proto::audit {
// A short query may end inside an observed k-mer. Virtual cursor nodes expose
// the remaining bases of THAT SAME observed k-mer before its ordinary edges.
// N padding is bookkeeping outside the emitted interval, never sequence data.
template<class Graph> class RootGraph {
    Graph& graph_;bool back_;
    struct Cursor {std::string label;Node next;char emitted;};
    std::map<Node,Cursor> cursors_;
public:
    RootGraph(Graph& graph,bool back):graph_(graph),back_(back){}
    size_t k()const{return graph_.k();}
    Path native_membership_roots(const Path& roots)const {
        Path result;
        // Reverse reachability cannot discover synthetic cursor predecessors
        // through native graph adjacency. Starting at the physical endpoint
        // instead is an optimistic membership relaxation, never an exclusion.
        for(auto node:roots) {
            auto it=cursors_.find(node);
            while(it!=cursors_.end()){node=it->second.next;it=cursors_.find(node);}
            result.push_back(node);
        }
        std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());
        return result;
    }
    Node add_root(Node node,const std::string& prefix) {
        if(prefix.empty())return node;
        auto word=graph_.label(node),expected=back_?word.substr(0,prefix.size()):word.substr(k()-prefix.size());
        if(back_)std::reverse(expected.begin(),expected.end());
        if(expected!=prefix)throw std::runtime_error("short-query cursor disagrees with observed k-mer");
        Node next=node;
        for(size_t remaining=1;remaining<=prefix.size();++remaining) {
            size_t consumed=prefix.size()-remaining;
            auto label=back_?word.substr(remaining)+std::string(remaining,'N'):std::string(remaining,'N')+word.substr(0,k()-remaining);
            Node id=missing-1-cursors_.size();cursors_.emplace(id,Cursor{label,next,prefix[consumed]});next=id;
        }
        return next;
    }
    Node find(const std::string& word)const{return graph_.find(word);}
    std::string label(Node v)const{auto it=cursors_.find(v);return it==cursors_.end()?graph_.label(v):it->second.label;}
    std::vector<Arc> neighbors(Node v,bool back)const {
        auto it=cursors_.find(v);if(it==cursors_.end())return graph_.neighbors(v,back);
        if(back!=back_)return {};return {{it->second.next,v,0}};
    }
    char extension_base(const Arc& arc,bool back)const {
        auto it=cursors_.find(arc.edge);if(it!=cursors_.end())return it->second.emitted;
        return graph_.extension_base(arc,back);
    }
    std::string spell(const Path& path)const {
        if(path.empty())return {};auto result=label(path.front());
        for(size_t i=1;i<path.size();++i) {
            auto from=path[i-1],to=path[i];bool valid=false;
            if(back_){auto it=cursors_.find(to);valid=it!=cursors_.end()&&it->second.next==from;}
            else {auto it=cursors_.find(from);valid=it!=cursors_.end()&&it->second.next==to;}
            if(!valid)for(auto arc:graph_.neighbors(from,false))if(arc.to==to)valid=true;
            if(!valid)throw std::runtime_error("disconnected cursor witness");
            result+=label(to).back();
        }
        return result;
    }
    double count(Node a,Node b)const {
        if(cursors_.count(a)||cursors_.count(b))return 0;return graph_.count(a,b);
    }
};
}
#endif
