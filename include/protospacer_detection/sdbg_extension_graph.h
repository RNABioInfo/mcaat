#ifndef MCAAT_SDBG_EXTENSION_GRAPH_H
#define MCAAT_SDBG_EXTENSION_GRAPH_H
#include "protospacer_detection/extension_engine.h"
#include "sdbg/sdbg.h"
#include <functional>

namespace mcaat::proto {
// Identity is a k-mer vertex (the last edge slot), while coverage belongs to the
// actual traversed (k+1)-mer edge. Never change global edge-validity flags.
class ExtensionGraph {
    const SDBG& graph_;
public:
    explicit ExtensionGraph(const SDBG& graph):graph_(graph) {}
    size_t k() const {return graph_.k();}
    std::string label(Node v) const {
        if(v>=graph_.size())throw std::runtime_error("vertex outside graph");
        std::vector<uint8_t> bases(k());graph_.GetLabel(v,bases.data());
        std::string s;for(auto b:bases){if(b<1 || b>4)throw std::runtime_error("non-DNA graph vertex");s+="ACGT"[b-1];}return s;
    }
    Node find(const std::string& s) const {
        if(s.size()!=k())return missing;
        std::vector<uint8_t> b;for(char c:s){auto i=std::string("ACGT").find(c);if(i==std::string::npos)return missing;b.push_back(i+1);}
        Node v=graph_.IndexBinarySearch(b.data());
        return v==missing || label(v)!=s?missing:v;
    }
    // BOSS vertices are ordered colexicographically. A suffix therefore forms
    // one contiguous range; do not enumerate all 4^(k-pattern_length) paddings.
    Path find_suffix(const std::string& suffix,size_t limit,bool& limited) const {
        if(suffix.empty() || suffix.size()>k() || suffix.find_first_not_of("ACGT")!=std::string::npos)return {};
        std::vector<uint8_t> bases(k());
        auto compare=[&](Node vertex) {
            graph_.GetLabel(vertex,bases.data());
            for(size_t j=0;j<suffix.size();++j) {
                int actual=bases[k()-1-j];int wanted=int(std::string("ACGT").find(suffix[suffix.size()-1-j]))+1;
                if(actual!=wanted)return actual<wanted?-1:1;
            }
            return 0;
        };
        Node lo=0,hi=graph_.size();
        while(lo<hi){Node mid=lo+(hi-lo)/2;if(compare(mid)<0)lo=mid+1;else hi=mid;}
        std::set<Node> unique;
        for(Node index=lo;index<graph_.size() && compare(index)==0;) {
            Node vertex=graph_.GetLastIndex(index);
            graph_.GetLabel(vertex,bases.data());std::string word;
            for(auto base:bases){if(base<1 || base>4){word.clear();break;}word+="ACGT"[base-1];}
            if(!word.empty()) {
                auto mapped=find(word);
                if(mapped!=missing && !unique.count(mapped)) {
                    if(unique.size()>=limit){limited=true;break;}
                    unique.insert(mapped);
                }
            }
            if(vertex<index)throw std::runtime_error("invalid BOSS vertex range");
            index=vertex+1;
        }
        return Path(unique.begin(),unique.end());
    }
    // Search observed BOSS labels with '.' wildcards and a substitution budget.
    // Work is bounded independently of the number of surviving graph paths.
    std::vector<std::pair<Node,size_t>> find_masked(const std::string& pattern,size_t budget,
            size_t max_work,size_t& work,bool& limited) const {
        if(pattern.size()!=k()||pattern.find_first_not_of("ACGT.")!=std::string::npos)
            throw std::runtime_error("invalid masked k-mer pattern");
        std::map<Node,size_t> result;std::vector<uint8_t> bases(k());
        auto spend=[&](){if(work>=max_work){limited=true;return false;}++work;return true;};
        std::function<void(Node,Node,size_t,size_t)> visit;
        visit=[&](Node lo,Node hi,size_t depth,size_t errors) {
            if(lo>=hi||!spend())return;
            if(depth==k()||hi-lo<=32) {
                for(Node index=lo;index<hi;) {
                    if(!spend())return;
                    Node v=graph_.GetLastIndex(index);graph_.GetLabel(v,bases.data());std::string word;
                    for(auto b:bases){if(b<1||b>4){word.clear();break;}word+="ACGT"[b-1];}
                    if(!word.empty()) {
                        size_t actual=0;for(size_t j=0;j<k();++j)actual+=pattern[j]!='.'&&pattern[j]!=word[j];
                        auto mapped=find(word);if(mapped!=missing&&actual<=budget)result.emplace(mapped,actual);
                    }
                    if(v<index)throw std::runtime_error("invalid BOSS range");index=v+1;
                }
                return;
            }
            const size_t at=k()-1-depth;
            auto bound=[&](unsigned base) {
                auto a=lo,b=hi;
                while(a<b){auto mid=a+(b-a)/2;graph_.GetLabel(mid,bases.data());if(bases[at]<base)a=mid+1;else b=mid;}
                return a;
            };
            for(unsigned base=1;base<=4&&!limited;++base) {
                size_t cost=errors+(pattern[at]!='.'&&pattern[at]!="ACGT"[base-1]);
                if(cost<=budget)visit(bound(base),bound(base+1),depth+1,cost);
            }
        };
        visit(0,graph_.size(),0,0);return {result.begin(),result.end()};
    }
    std::vector<Arc> neighbors(Node v,bool backward) const {
        std::vector<Arc> result;if(v>=graph_.size())return result;
        auto add=[&](Node e,Node to) {
            auto b=graph_.GetW(e);if(b>4)b-=4;
            if(b>=1 && b<=4 && graph_.IsValidEdge(e))result.push_back({to,e,double(graph_.EdgeMultiplicity(e))});
        };
        if(!backward) {
            Node e=v;do {
                auto w=graph_.GetW(e);if(w>4)w-=4;
                if(w>=1 && w<=4 && graph_.IsValidEdge(e))add(e,graph_.Forward(e));
                --e;
            }while(e!=missing && !graph_.IsLastOrTip(e));
        } else {
            // Like ComputeIncomings, but a sink's invalid sentinel edge does not
            // invalidate the real incoming edges of its k-mer vertex.
            Node first=graph_.Backward(v);if(first>=graph_.size())return result;
            auto c=graph_.GetW(first);unsigned last=graph_.IsLastOrTip(first);
            add(first,graph_.GetLastIndex(first));
            for(Node e=first+1;last<5 && e<graph_.size();++e) {
                last+=graph_.IsLastOrTip(e);auto w=graph_.GetW(e);
                if(w==c)break;
                if(w==c+4)add(e,graph_.GetLastIndex(e));
            }
        }
        std::sort(result.begin(),result.end(),[](auto a,auto b){return a.to<b.to;});
        return result;
    }
    Node edge(Node a,Node b) const {for(auto e:neighbors(a,false))if(e.to==b)return e.edge;return missing;}
    char extension_base(const Arc& arc,bool backward) const {
        if(backward)return label(arc.to).front();
        auto b=graph_.GetW(arc.edge);if(b>4)b-=4;
        if(b<1 || b>4)throw std::runtime_error("non-DNA extension edge");
        return "ACGT"[b-1];
    }
    double count(Node a,Node b) const {for(auto e:neighbors(a,false))if(e.to==b)return e.count;throw std::runtime_error("missing path edge");}
    Path trace(const std::string& s) const {
        if(s.size()<k())return {};
        Node v=find(s.substr(0,k()));if(v==missing)return {};
        Path path{v};
        for(size_t i=k();i<s.size();++i) {
            Node next=missing;
            for(auto e:neighbors(v,false)) {
                auto code=graph_.GetW(e.edge);if(code>4)code-=4;
                if(code>=1 && code<=4 && "ACGT"[code-1]==s[i]){next=e.to;break;}
            }
            if(next==missing)return {};
            path.push_back(next);v=next;
        }
        if(spell(path)!=s)throw std::runtime_error("mapped path does not spell input");
        return path;
    }
    std::string spell(const Path& path) const {
        if(path.empty())return {};
        std::string s=label(path.front());
        for(size_t i=1;i<path.size();++i) {
            Node e=edge(path[i-1],path[i]);if(e==missing)throw std::runtime_error("disconnected extension path");
            auto b=graph_.GetW(e);if(b>4)b-=4;
            s+="ACGT"[b-1];
        }
        return s;
    }
    std::pair<std::set<Node>,std::set<Node>> array_context(const std::string& repeat,const std::string& spacer) const {
        std::set<Node> edges,pure_repeat;
        for(const auto& pair:std::vector<std::pair<std::string,std::string>>{{repeat,spacer},{reverse_complement(repeat),reverse_complement(spacer)}}) {
            auto seq=pair.first+pair.second+pair.first;
            size_t begin=pair.first.size(),end=begin+pair.second.size();
            Node previous=missing;
            for(size_t i=0;i+k()<=seq.size();++i) {
                Node v=missing;
                if(previous!=missing)for(auto arc:neighbors(previous,false)) {
                    auto b=graph_.GetW(arc.edge);if(b>4)b-=4;
                    if(b>=1 && b<=4 && "ACGT"[b-1]==seq[i+k()-1]){v=arc.to;break;}
                }
                if(v==missing)v=find(seq.substr(i,k()));
                if(i && previous!=missing && v!=missing && !(i-1>=begin && i+k()<=end)) {
                    Node e=edge(previous,v);if(e!=missing) {
                        edges.insert(e);
                        if(i+k()<=begin || i-1>=end)pure_repeat.insert(e);
                    }
                }
                previous=v;
            }
        }
        return {edges,pure_repeat};
    }
};
int run_extension_module(int argc,char** argv);
} // namespace mcaat::proto
#endif
