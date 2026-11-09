#include "protospacer_detection/extension_engine.h"
#include "protospacer_detection/extension_evidence.h"
#include <iostream>
#include <random>
using namespace mcaat::proto;
size_t checks=0;
void check(bool condition,const std::string& why){++checks;if(!condition)throw std::runtime_error(why);}
template<class F> void invalid(F f){bool failed=false;try{f();}catch(const std::exception&){failed=true;}check(failed,"invalid input accepted");}
struct Graph {
    size_t k_=11;std::map<std::string,Node> ids;std::vector<std::string> labels;
    std::map<std::pair<Node,Node>,Arc> arcs;
    size_t k()const{return k_;}
    Node node(const std::string& s){auto [it,added]=ids.emplace(s,ids.size());if(added)labels.push_back(s);return it->second;}
    Path add(const std::string& sequence,double count=12) {
        Path p;for(size_t i=0;i+k()<=sequence.size();++i)p.push_back(node(sequence.substr(i,k())));
        for(size_t i=1;i<p.size();++i){auto [it,added]=arcs.emplace(std::make_pair(p[i-1],p[i]),Arc{p[i],arcs.size()+1,0});it->second.count+=count;}
        return p;
    }
    std::vector<Arc> neighbors(Node v,bool backward)const {
        std::vector<Arc> result;for(auto [ends,a]:arcs)if((backward?ends.second:ends.first)==v){a.to=backward?ends.first:ends.second;result.push_back(a);}return result;
    }
    std::string label(Node v)const{return labels.at(v);}
    double count(Node a,Node b)const{return arcs.at({a,b}).count;}
    std::string spell(const Path& p)const {
        if(p.empty())return {};std::string s=label(p[0]);
        for(size_t i=1;i<p.size();++i){if(!arcs.count({p[i-1],p[i]}))throw std::runtime_error("disconnected");s+=label(p[i]).back();}return s;
    }
};
std::mt19937 random_dna(8237);
std::string dna(size_t n){std::string s;while(n--)s+="ACGT"[random_dna()%4];return s;}
struct Fixture {
    Graph g;Options o;std::string repeat=dna(24),spacer=dna(32),left=dna(100),right=dna(150);
    Path seed;std::set<Node> masked;
    Fixture(bool both=true) {
        auto array=g.add(repeat+spacer+repeat,30);
        for(size_t i=1;i<array.size();++i)if(i+g.k()<=repeat.size() || i-1>=repeat.size()+spacer.size())masked.insert(g.arcs.at({array[i-1],array[i]}).edge);
        g.add((both?left:std::string())+spacer+right,12);
        for(size_t i=0;i+g.k()<=spacer.size();++i)seed.push_back(g.ids.at(spacer.substr(i,g.k())));
        o.max_bp=500;o.coverage_window=3;o.max_coverage_penalty=20;
    }
    Seed full() {auto found=isolate(g,seed,masked,o);check(!found.seeds.empty(),"seed detected");check(found.seeds[0].begin==0 && found.seeds[0].end==spacer.size(),"full match ranked before fragments");return found.seeds[0];}
};
void test_extension() {
    Fixture f;auto seed=f.full();auto found=reconstruct(f.g,seed,f.masked,f.o);
    check(std::any_of(found.begin(),found.end(),[&](const auto& c){return c.sequence==f.left+f.spacer+f.right;}),"exact bidirectional reconstruction");
    for(const auto& c:found) {
        check(c.sequence.substr(c.seed_begin,c.seed_end-c.seed_begin)==f.spacer,"seed coordinates spell spacer");
        check(c.sequence.size()==c.path.size()+f.g.k()-1,"path length accounting");
        if(c.sequence==f.left+f.spacer+f.right)check(c.left.coverage.mean()==12 && c.right.coverage.mean()==12,"shared spacer depth excluded from flank baseline");
        check(c.sequence.size()<=f.o.max_bp,"length budget respected");
    }
    Fixture one(false);auto partial=reconstruct(one.g,one.full(),one.masked,one.o);
    check(std::any_of(partial.begin(),partial.end(),[&](auto& c){return c.sequence.size()>=one.spacer.size()+one.right.size() && c.sequence.substr(c.seed_begin)==one.spacer+one.right && c.right_length>0;}),"one-sided sequence retained, including shared array-boundary prefix");
    f.o.max_bp=80;found=reconstruct(f.g,seed,f.masked,f.o);
    check(std::all_of(found.begin(),found.end(),[](auto& c){return c.sequence.size()<=80;}),"total, not per-side, bp budget");
    f.o.max_states=2;found=reconstruct(f.g,seed,f.masked,f.o);
    check(!found.empty() && found[0].search_limited && found[0].expanded<=4,"work limit preserves seed and reports truncation");
}
void test_branches() {
    Fixture f;auto wrong=dna(70);f.g.add(f.spacer+f.right.substr(0,70)+wrong,90);
    // Reset the shared correct flank to 12 so the branch test isolates a depth jump.
    auto correct=f.spacer+f.right;
    for(size_t i=1;i+f.g.k()<=correct.size();++i) {
        auto a=f.g.ids.at(correct.substr(i-1,f.g.k())),b=f.g.ids.at(correct.substr(i,f.g.k()));f.g.arcs.at({a,b}).count=12;
    }
    auto result=reconstruct(f.g,f.full(),f.masked,f.o);
    check(std::any_of(result.begin(),result.end(),[&](auto& c){return c.sequence==f.left+f.spacer+f.right;}),"coverage-compatible branch survives");
    check(std::none_of(result.begin(),result.end(),[&](auto& c){return c.sequence.find(wrong)!=std::string::npos;}),"high-depth branch is stopped rather than rewarded");
    auto found=isolate(f.g,f.seed,f.masked,f.o);
    f.o.max_seeds=1;auto capped=isolate(f.g,f.seed,f.masked,f.o);
    check(capped.seeds.size()==1 && capped.eligible==found.eligible,"seed cap reports total eligible intervals");
    // A second target diverging inside the spacer gives overlapping intervals;
    // neither may be removed because it shares seed nodes with another.
    f.g.add(dna(30)+f.spacer.substr(4,22)+dna(30),8);
    found=isolate(f.g,f.seed,f.masked,Options{});
    check(found.seeds.size()>1,"overlapping seed intervals retained");
}
void test_coverage() {
    check(Coverage::loss(12,12,0)==0,"Poisson saturated mean");
    check(Coverage::loss(13,12,0)<Coverage::loss(90,12,0),"compatibility, not maximum count");
    check(Coverage::loss(90,12,.5)<Coverage::loss(90,12,0),"overdispersion widens count tolerance");
    Options o;o.coverage_window=3;Coverage c;
    for(double x:{2,20,60})c.observe(x,o);
    check(c.dispersion()>0,"local dispersion estimated");
    for(double x:{0.,1.,3.,10.,50.,100.})for(double m:{1.,10.,50.})for(double d:{0.,.1,1.})
        check(std::isfinite(Coverage::loss(x,m,d)) && Coverage::loss(x,m,d)>=0,"finite nonnegative coverage loss");
}
void test_join_and_repeated_seed() {
    Graph g;Options o;o.max_bp=500;
    auto spacer=dna(32),shared=dna(45);
    auto cycle=g.add(shared+spacer+shared);
    Path path(cycle.begin()+shared.size(),cycle.begin()+shared.size()+spacer.size()-g.k()+1);
    auto seeds=isolate(g,path,{},o);check(!seeds.seeds.empty(),"cycle seed found");
    auto joined=reconstruct(g,seeds.seeds.front(),{},o);
    check(joined.size()==1 && joined[0].sequence==spacer,"intersecting flanks cannot create a duplicated region");
    check(joined[0].left.stop=="incompatible_flanks" && joined[0].expanded>0,"failed join retains seed and actual work");
    Graph repeated;auto motif=dna(18);spacer=dna(16)+motif+dna(13)+motif+dna(16);
    auto left=dna(40),right=dna(50);auto full=repeated.add(left+spacer+right);
    path=Path(full.begin()+left.size(),full.begin()+left.size()+spacer.size()-repeated.k()+1);
    check(std::set<Node>(path.begin(),path.end()).size()<path.size(),"fixture contains repeated seed k-mers");
    seeds=isolate(repeated,path,{},o);check(!seeds.seeds.empty(),"repeated-kmer seed retained");
    joined=reconstruct(repeated,seeds.seeds.front(),{},o);
    check(std::any_of(joined.begin(),joined.end(),[&](const auto& c){return c.sequence==left+spacer+right;}),"repeated seed does not invalidate disjoint target flanks");
}
void test_reads_and_pam() {
    ReadEvidence evidence(4,1000);
    evidence.index({{"AACCGGTTAATGCC",{"CCGGTT"}},{"AACCGGTTCCAGCC",{"CCGGTT"}}});
    evidence.read("AACCGGTTA");evidence.read(reverse_complement("AACCGGTTA"));
    check(evidence.junction_count("CCGGTT")==2,"read junction support is strand invariant");
    evidence.pair("AACCGG",reverse_complement("AATGCC"),0,100);
    check(evidence.regions[0].unique_pairs==1 && evidence.regions[1].pairs==0,"full mate sequences choose supported alternative");
    evidence.pair("AACCGG",reverse_complement("CCGGTT"),0,100);
    check(evidence.regions[0].pairs==2 && evidence.regions[1].pairs==1 && evidence.regions[1].unique_pairs==0,"ambiguous pairs not counted as unique evidence");
    evidence.pair("AACCGG","AATGCC",0,100);
    check(evidence.regions[0].pairs==2,"wrong mate orientation rejected");
    evidence.pair("AACCGG",reverse_complement("AATGCC"),15,100);
    check(evidence.regions[0].pairs==2,"insert size bounds enforced");
    invalid([]{ReadEvidence r(3,1);r.index({{"ACGTACGT",{}}});});
    invalid([]{ReadEvidence r(3,1000,1);r.index({{"ACGTACGTACGTACGT",{}}});r.pair("ACGT","ACGT",0,100);});
    PamRule rule{"Array_1","I-E","+","3prime","NGG"};
    check(check_pam("AACGTAGG",1,5,false,rule,true)=="match","PAM matches on specified side");
    check(check_pam("AACGTAGG",1,5,false,rule,false)=="unknown_partial_match","partial spacer cannot imply PAM position");
    check(check_pam("AACGTA",1,5,false,rule,true)=="unknown_short_flank","missing flank is unknown");
    rule.strand="?";check(check_pam("AACGTAGG",1,5,false,rule,true)=="unknown_orientation","orientation is not guessed");
    rule.strand="-";check(check_pam(reverse_complement("AACGTAGG"),3,7,false,rule,true)=="match","reverse orientation PAM coordinates");
    rule.side="none";rule.motif="-";check(check_pam("AACGTA",1,5,false,rule,false)=="not_applicable","PAM-inapplicable model is explicit");
    for(char c:std::string("ACGTRYSWKMBDHVN"))check(iupac(c,'A') || iupac(c,'C') || iupac(c,'G') || iupac(c,'T'),"IUPAC motif supported");
    invalid([]{iupac('Z','A');});
}
int main() {
    try {
        test_extension();test_branches();test_coverage();test_reads_and_pam();test_join_and_repeated_seed();
        for(auto model:{"auto","poisson","negative-binomial"}){Options o;o.coverage_model=model;o.validate();}
        invalid([]{Options o;o.beam_width=0;o.validate();});
        invalid([]{Options o;o.max_coverage_penalty=std::nan("");o.validate();});
        invalid([]{Options o;o.coverage_model="bad";o.validate();});
        std::cout<<checks<<" protospacer checks passed\n";
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
