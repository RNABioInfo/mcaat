#include "../bench/reference_graph_oracle.h"
#include <iostream>
#include <random>
using namespace mcaat::proto;
struct Graph {
    size_t k_=23;
    size_t k()const{return k_;}
    std::vector<std::string> labels;std::map<std::string,Node> ids;std::map<std::pair<Node,Node>,Arc> arcs;
    Path add(const std::string& sequence){Path p;for(size_t i=0;i+k()<=sequence.size();++i){auto s=sequence.substr(i,k());auto [it,ok]=ids.emplace(s,ids.size());if(ok)labels.push_back(s);p.push_back(it->second);if(p.size()>1)arcs[{p[p.size()-2],p.back()}]={p.back(),arcs.size(),12};}return p;}
    std::string label(Node v)const{return labels.at(v);}
    std::vector<Arc> neighbors(Node v,bool back)const{std::vector<Arc> r;for(auto [e,a]:arcs)if((back?e.second:e.first)==v){a.to=back?e.first:e.second;r.push_back(a);}return r;}
    char extension_base(Arc a,bool back)const{return back?label(a.to).front():label(a.to).back();}
    std::string spell(const Path& p)const{auto s=label(p.front());for(size_t i=1;i<p.size();++i){if(!arcs.count({p[i-1],p[i]}))throw std::runtime_error("disconnected");s+=label(p[i]).back();}return s;}
};
// Independent string DP over every prefix of a single nonbranching path.
size_t optimal(const std::string& ref,const std::string& seq){std::vector<size_t> row(seq.size()+1),next(row.size());for(size_t j=0;j<row.size();++j)row[j]=j*9;for(size_t i=1;i<=ref.size();++i){next[0]=i*10;for(size_t j=1;j<row.size();++j)next[j]=std::min({row[j]+10,next[j-1]+9,row[j-1]+(ref[i-1]==seq[j-1]?0:10)});row.swap(next);}return *std::min_element(row.begin(),row.end());}
int main(){try{
    std::mt19937 rng(111129);auto dna=[&](size_t n){std::string s;while(n--)s+="ACGT"[rng()%4];return s;};
    size_t cases=0;
    for(int variant=0;variant<24;++variant){auto seed=dna(23),ref=dna(100),seq=ref;
        if(variant<16)for(int m=0;m<variant;++m){size_t p=m*6;seq[p]="ACGT"[(std::string("ACGT").find(seq[p])+1)%4];}
        else if(variant<20)seq.insert(40,dna((variant-15)*3));else seq.erase(35,(variant-19)*3);
        Graph g;auto p=g.add(seed+seq);auto result=audit::align(g,p.front(),ref,false,{100},200000);
        bool expected=optimal(ref,seq)<=100;
        if(result.found.count(100)!=expected || result.limited)throw std::runtime_error("oracle disagrees with independent string DP: "+std::to_string(variant));
        if(expected){auto w=result.found.at(100);if(w.matches*100<90*w.columns || seq.substr(0,w.sequence.size())!=w.sequence)throw std::runtime_error("bad witness");}
        auto backward_ref=ref;std::reverse(backward_ref.begin(),backward_ref.end());
        Graph bg;auto bp=bg.add(seq+seed);auto back=audit::align(bg,bp.back(),backward_ref,true,{100},200000);
        std::reverse(seq.begin(),seq.end());if(back.found.count(100)!=(optimal(backward_ref,seq)<=100))throw std::runtime_error("backward oracle mismatch");
        std::reverse(seq.begin(),seq.end());
        auto bounded=audit::align(g,p.front(),ref,false,{24,100},120);
        if(!bounded.found.count(100) && bounded.lowest_unexplored_cost>100 && expected)
            throw std::runtime_error("false certificate of absent qualifying path");
        if(!bounded.found.count(24) && bounded.lowest_unexplored_cost>24 && optimal(ref.substr(0,24),seq)<=24)
            throw std::runtime_error("false short-prefix absence certificate");
        ++cases;
    }
    Graph g;auto p=g.add(dna(160));auto small=audit::align(g,p.front(),dna(100),false,{100},2);
    if(!small.limited||small.exhausted)throw std::runtime_error("capped search mislabeled as absence");
    // Exhaustively enumerate short walks in branching/merging/cyclic graphs.
    // This independently checks product-state dominance, not just linear paths.
    for(size_t variant=0;variant<12;++variant) {
        Graph graph;graph.k_=3;auto native=dna(14);auto path=graph.add(native);
        for(size_t j=0;j<4;++j)graph.add(dna(10));
        auto ref=variant%2?native.substr(3,10):dna(10);
        bool expected=false;
        std::function<void(Node,std::string)> visit=[&](Node v,std::string sequence) {
            if(optimal(ref,sequence)<=ref.size()){expected=true;return;}
            if(sequence.size()>=11 || expected)return;
            for(auto arc:graph.neighbors(v,false))visit(arc.to,sequence+graph.extension_base(arc,false));
        };
        visit(path.front(),"");
        auto actual=audit::align(graph,path.front(),ref,false,{10},50000);
        if(actual.limited || bool(actual.found.count(10))!=expected)
            throw std::runtime_error("oracle disagrees with exhaustive branching walk enumeration");
    }
    std::cout<<cases<<" forward/backward DP cases and 12 exhaustive branching graphs passed; capped search stays unresolved\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
