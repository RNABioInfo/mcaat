#include "protospacer_detection/retained_target.h"
#include <fstream>
#include <iostream>
#include <sstream>
using namespace mcaat::proto;
// Exact observed SDBG edges and counts, restricted to diagnostic regions.
struct Graph {
    std::map<std::string,Node> ids;std::vector<std::string> labels;
    std::map<Node,std::vector<Arc>> out,in;
    size_t k()const{return 23;}
    Node node(const std::string& s){auto [it,added]=ids.emplace(s,ids.size());if(added)labels.push_back(s);return it->second;}
    void add(const std::string& a,const std::string& b,double count) {
        if(a.size()!=k() || b.size()!=k() || a.substr(1)!=b.substr(0,k()-1))throw std::runtime_error("invalid fixture edge");
        auto from=node(a),to=node(b);auto id=from*4+std::string("ACGT").find(b.back());
        out[from].push_back({to,id,count});in[to].push_back({from,id,count});
    }
    std::vector<Arc> neighbors(Node v,bool back)const{auto& m=back?in:out;auto it=m.find(v);return it==m.end()?std::vector<Arc>{}:it->second;}
    std::string label(Node v)const{return labels.at(v);}
    char extension_base(const Arc& a,bool back)const{return back?label(a.to).front():label(a.to).back();}
    Node edge(Node a,Node b)const{for(auto arc:neighbors(a,false))if(arc.to==b)return arc.edge;return missing;}
    Path trace(const std::string& seq)const {
        Path p;for(size_t i=0;i+k()<=seq.size();++i) {
            auto it=ids.find(seq.substr(i,k()));if(it==ids.end())return {};
            if(!p.empty() && edge(p.back(),it->second)==missing)return {};
            p.push_back(it->second);
        }return p;
    }
    std::string spell(const Path& p)const {
        if(p.empty())return {};auto s=label(p.front());
        for(size_t i=1;i<p.size();++i){if(edge(p[i-1],p[i])==missing)throw std::runtime_error("disconnected");s+=label(p[i]).back();}return s;
    }
    std::pair<std::set<Node>,std::set<Node>> array_context(const std::string& repeat,const std::string& spacer)const {
        std::set<Node> edges;auto seq=repeat+spacer+repeat;
        for(size_t i=0;i+k()<seq.size();++i) {
            if(i>=repeat.size() && i+k()+1<=repeat.size()+spacer.size())continue;
            auto p=trace(seq.substr(i,k()+1));if(!p.empty())edges.insert(edge(p[0],p[1]));
        }return {edges,{}};
    }
};
int main(int argc,char** argv) {
    try {
        if(argc!=2)throw std::runtime_error("fixture directory required");
        std::ifstream edges(std::string(argv[1])+"/edges.tsv"),cases(std::string(argv[1])+"/cases.tsv");
        if(!edges || !cases)throw std::runtime_error("missing real graph fixtures");
        Graph graph;std::string line;std::getline(edges,line);
        while(std::getline(edges,line)){std::istringstream row(line);std::string a,b;double count;row>>a>>b>>count;graph.add(a,b,count);}
        std::getline(cases,line);size_t tested=0;
        while(std::getline(cases,line)) {
            std::istringstream row(line);std::string id,seed,spacer,repeat,expected;size_t begin,end;
            row>>id>>seed>>spacer>>repeat>>begin>>end>>expected;
            auto result=extend_retained_target(graph,seed,{{spacer,repeat,begin,end}},Options{});
            if(expected=="-") {
                if(result.has_extension)throw std::runtime_error(id+": escaped through source repeat");
            }else if(!result.has_extension || result.selected.sequence!=expected) {
                throw std::runtime_error(id+": failed reference-supported continuation; observed length="+std::to_string(result.selected.sequence.size()));
            }
            ++tested;
        }
        if(tested!=2)throw std::runtime_error("missing real graph cases");
        std::cout<<tested<<" observed graph cases passed\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
