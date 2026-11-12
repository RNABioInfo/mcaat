#include "../bench/complete_target_placements.h"
#include <random>
#include <iostream>
using namespace mcaat::proto;
struct Graph {
    size_t k_=23;
    size_t k()const{return k_;}
    std::vector<std::string> labels;std::map<std::string,Node> ids;std::map<std::pair<Node,Node>,Arc> arcs;
    Path add(const std::string& sequence){Path p;for(size_t i=0;i+k()<=sequence.size();++i){auto s=sequence.substr(i,k());auto [it,ok]=ids.emplace(s,ids.size());if(ok)labels.push_back(s);p.push_back(it->second);if(p.size()>1)arcs.emplace(std::make_pair(p[p.size()-2],p.back()),Arc{p.back(),arcs.size(),12});}return p;}
    std::string label(Node v)const{return labels.at(v);}
    Node find(const std::string& s)const{auto it=ids.find(s);return it==ids.end()?missing:it->second;}
    std::vector<Arc> neighbors(Node v,bool back)const{std::vector<Arc> r;for(auto [e,a]:arcs)if((back?e.second:e.first)==v){a.to=back?e.first:e.second;r.push_back(a);}return r;}
    char extension_base(Arc a,bool back)const{return back?label(a.to).front():label(a.to).back();}
    std::string spell(const Path& p)const{auto s=label(p.front());for(size_t i=1;i<p.size();++i){if(!arcs.count({p[i-1],p[i]}))throw std::runtime_error("disconnected");s+=label(p[i]).back();}return s;}
    double count(Node,Node)const{return 12;}
    std::vector<std::pair<Node,size_t>> find_masked(const std::string& pattern,size_t budget)const {
        std::vector<std::pair<Node,size_t>> result;
        for(size_t v=0;v<labels.size();++v){size_t errors=0;for(size_t i=0;i<k();++i)errors+=pattern[i]!='.'&&pattern[i]!=labels[v][i];if(errors<=budget)result.emplace_back(v,errors);}
        return result;
    }
};

bool accepts(const std::string& query,const std::string& observed) {
    for(size_t b=0;b<query.size();++b)for(size_t e=b+1;e<=query.size();++e) {
        if((e-b)*100<query.size()*85)continue;
        size_t same=0;for(size_t i=b;i<e;++i)same+=query[i]==observed[i];
        if(same*100>=(e-b)*90)return true;
    }
    return false;
}
int main(){try {
    std::mt19937 rng(319871);auto dna=[&](size_t n){std::string s;while(n--)s+="ACGT"[rng()%4];return s;};
    size_t checked=0;
    for(size_t n:{15,17,19,22,23,26,30,35,38,42})for(size_t trial=0;trial<5;++trial) {
        Graph graph;auto sequence=dna(90);graph.add(sequence);
        auto variant=sequence;for(size_t i=10;i<70;i+=11)variant[i]="ACGT"[(std::string("ACGT").find(variant[i])+1)%4];graph.add(variant);
        auto query=sequence.substr(20,n);if(trial==1)query[0]='A';if(trial==2)query[n/2]='C';
        if(trial==3)query=dna(n);if(trial==4)query=variant.substr(18,n);
        std::set<std::pair<Node,std::string>> expected;
        for(Node v=0;v<graph.labels.size();++v) {
            if(n<graph.k()) {
                auto label=graph.label(v);
                for(size_t offset=0;offset+n<=label.size();++offset)if(accepts(query,label.substr(offset,n)))expected.emplace(v,label.substr(offset+n));
            }else {
                std::function<void(Path)> visit=[&](Path nodes){auto spelled=graph.spell(nodes);
                    if(spelled.size()==n){if(accepts(query,spelled))expected.emplace(nodes.back(),"");return;}
                    for(auto arc:graph.neighbors(nodes.back(),false)){auto next=nodes;next.push_back(arc.to);visit(next);}
                };visit(Path{v});
            }
        }
        auto result=audit::complete_target_placements(graph,query);std::set<std::pair<Node,std::string>> observed;
        for(auto [key,p]:result.endpoints){observed.insert(key);if(!accepts(query,p.observed))throw std::runtime_error("invalid witness");}
        if(expected!=observed)throw std::runtime_error("complete placement mismatch n="+std::to_string(n)+" trial="+std::to_string(trial));
        ++checked;
    }
    std::cout<<checked<<" complete placement comparisons against exhaustive graph walks passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
