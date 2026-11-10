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
    char extension_base(const Arc& a,bool backward)const{return backward?label(a.to).front():label(a.to).back();}
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
void test_array_context() {
    const std::string repeat="CGGTTCATCCCCGCTTGCGCGGGGAACTC";
    const std::string left="CGGTTCATCCCCGCTGGCGCGGGGAACTC";
    const std::string right="CGGGTAATCCCCGCTTGCGCGGGGAACTC";
    check(context_state(left,{{"",repeat}},.85)==ContextState::array,"observed repeat with one substitution rejected");
    check(context_state(right,{{"",repeat}},.85)==ContextState::array,"observed repeat with two substitutions rejected");
    check(context_state(repeat.substr(0,10)+"A"+repeat.substr(10),{{"",repeat}},.85)==ContextState::array,"repeat insertion handled");
    check(context_state(repeat.substr(0,10)+repeat.substr(11),{{"",repeat}},.85)==ContextState::array,"repeat deletion handled");
    check(context_state(repeat.substr(0,12),{{"",repeat}},.85)==ContextState::pending,"short repeat-compatible flank unresolved");
    check(context_state("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",{{"",repeat}},.85)==ContextState::external,"unrelated target flank retained");
    check(context_state("ACGT"+left,{{"ACGT",repeat}},.85)==ContextState::array,"partial seed continuing source spacer reaches its repeat");
    check(context_state("ATGT"+left,{{"ACGT",repeat}},.85)==ContextState::array,"related spacer variant remains array context");
    check(context_state("ATGT"+std::string(30,'A'),{{"ACGT",repeat}},.85)==ContextState::external,"partial target with unrelated flank retained");
    check(context_state(left,{{"",repeat}},1)==ContextState::external,"repeat similarity is user adjustable");
    auto backwards=flank_contexts("ACGT",2,4,{repeat},true);
    auto reversed=left;std::reverse(reversed.begin(),reversed.end());
    check(context_state("CA"+reversed,backwards,.85)==ContextState::array,"backward context reverses bases without complementing twice");
    Graph graph;graph.k_=23;Options o;o.max_bp=500;
    auto spacer=std::string("ACTACTTCGCTGGAAACGCTAGAGGTGCGTTTT");
    auto host=graph.add(dna(70)+left+spacer+right+dna(70));
    Path path;for(size_t i=0;i+graph.k()<=spacer.size();++i)path.push_back(graph.ids.at(spacer.substr(i,graph.k())));
    auto seeds=isolate(graph,path,{},o);check(!seeds.seeds.empty(),"variant host seed present");
    auto lc=flank_contexts(spacer,0,spacer.size(),{repeat},true),rc=flank_contexts(spacer,0,spacer.size(),{repeat},false);
    auto candidates=reconstruct(graph,seeds.seeds[0],{},o,nullptr,&lc,&rc);
    check(candidates.size()==1 && candidates[0].sequence==spacer,"source array cannot become an extended target");
    check(candidates[0].left.context=="array" && candidates[0].right.context=="array","both variant repeats recognized");
    auto zero_budget=o;zero_budget.max_bp=spacer.size();
    auto untested=reconstruct(graph,seeds.seeds[0],{},zero_budget,nullptr,&lc,&rc);
    check(untested[0].left.context=="incomplete" || untested[0].right.context=="incomplete","zero extension budget is unresolved, not a negative result");
    auto target_left=dna(60),target_right=dna(80);graph.add(target_left+spacer+target_right);
    seeds=isolate(graph,path,{},o);candidates=reconstruct(graph,seeds.seeds[0],{},o,nullptr,&lc,&rc);
    check(std::any_of(candidates.begin(),candidates.end(),[&](const auto& c){return c.sequence==target_left+spacer+target_right;}),"true target for the same spacer survives source-array filtering");
    // An exact prefix longer than k can be shared by a genuine flank. Only the
    // complete repeat similarity decision may block it, not a global edge mask.
    Graph shared;shared.k_=11;auto r=dna(40),s=dna(32),suffix=r.substr(0,18)+dna(50);
    auto true_path=shared.add(dna(50)+s+suffix);shared.add(r+s+r);
    path.clear();for(size_t i=0;i+shared.k()<=s.size();++i)path.push_back(shared.ids.at(s.substr(i,shared.k())));
    seeds=isolate(shared,path,{},o);lc=flank_contexts(s,0,s.size(),{r},true);rc=flank_contexts(s,0,s.size(),{r},false);
    candidates=reconstruct(shared,seeds.seeds[0],{},o,nullptr,&lc,&rc);
    check(std::any_of(candidates.begin(),candidates.end(),[&](const auto& c){return c.sequence==shared.spell(true_path);}),"long shared repeat prefix can diverge to a real target");
}
void test_interval_compaction() {
    RetainedIntervalContexts index;
    check(index.replacement("AACCGGTT",2,0,4,"full").empty(),"maximal interval retained");
    check(index.replacement("AACCGGTT",3,1,3,"partial")=="full","nested interval with identical context consolidated");
    check(index.replacement("AACCGGTA",3,1,3,"other").empty(),"different target path retained");
    check(index.replacement("AACCGGTT",4,1,3,"shifted").empty(),"different alignment retained");
    check(index.replacement("AACCGGTT",5,3,6,"overlap").empty(),"overlap without containment retained");
}
void test_stability_selection() {
    Options o;o.max_coverage_penalty=100000;
    auto profile=[&](std::vector<double> values) {Coverage c;for(double v:values)c.observe(v,o);return c;};
    auto flat=profile(std::vector<double>(20,20)),high=profile(std::vector<double>(20,200));
    auto varied=profile({5,40,5,40,5,40,5,40,5,40,5,40,5,40,5,40,5,40,5,40});
    check(stability(flat,flat).score()>stability(varied,varied).score(),"stable multiplicity outranks variable multiplicity");
    check(stability(flat,flat).score()>stability(flat,high).score(),"individually stable but incompatible flank levels penalized");
    check(std::abs(stability(flat,flat).score()-stability(high,high).score())<1e-12,"higher absolute multiplicity is not itself rewarded");
    auto short_flat=profile({20}),empty=profile({});
    check(stability(flat,flat).score()>stability(short_flat,short_flat).score(),"one observation cannot masquerade as strongest stable evidence");
    check(stability(flat,flat).score()>stability(empty,empty).score(),"absent coverage cannot earn perfect stability");
    auto blocks=profile({10,10,10,10,30,30,30,30}),alternating=profile({10,30,10,30,10,30,10,30});
    check(std::abs(blocks.profile.variance()-alternating.profile.variance())<1e-12,"equal count distributions have equal variance");
    check(stability(blocks,blocks).score()>stability(alternating,alternating).score(),"repeated jumps penalized beyond overall variance");
    auto noisy_baseline=profile({1,20,50});
    check(noisy_baseline.profile.variance()>0 && noisy_baseline.scored==0,"initial baseline variation participates in final ranking");
    Candidate a,b;a.sequence="ACGT";b.sequence="ACGA";a.left.coverage=a.right.coverage=flat;b.left.coverage=flat;b.right.coverage=high;
    a.score=stability(a.left.coverage,a.right.coverage).score();b.score=stability(b.left.coverage,b.right.coverage).score();
    check(better_prediction(a,30,b,30,o),"cross-flank stability drives prediction comparison");
    b.pair_support=2;check(better_prediction(b,30,a,30,o),"supplied physical linkage remains explicit ranking evidence");
    Graph g;g.k_=23;auto spacer=dna(34),l1=dna(260),r1=dna(280),l2=dna(270),r2=dna(250);
    g.add(l1+spacer+r1,12);g.add(l2+spacer+r2,80);
    Path seed;for(size_t i=0;i+g.k()<=spacer.size();++i)seed.push_back(g.ids.at(spacer.substr(i,g.k())));
    o.max_bp=1000;auto seeds=isolate(g,seed,{},o);auto candidates=reconstruct(g,seeds.seeds.front(),{},o);
    check(candidates.size()>=4,"fixture explores both matching and crossed flank combinations");
    auto best=std::min_element(candidates.begin(),candidates.end(),[&](const auto& x,const auto& y){return better_prediction(x,spacer.size(),y,spacer.size(),o);});
    check(best!=candidates.end() && (best->sequence==l1+spacer+r1 || best->sequence==l2+spacer+r2),"stable matching flank pair wins over chimeric low/high combinations");
}
void test_adapter_branches() {
    const std::string adapter="AGATCGGAAGAGC";
    const std::string observed_seed="AAACCATCGTTAGTTGTGATGCTACGAAT";
    // The failing exported path changes from a Shah-reference-matched seed
    // directly into this adapter prefix. Synthetic flanks isolate the mechanism.
    const std::string observed_tail="AGATCGGAAGAGCACAAGTCTGAACACCAGTCACAGATCAGTATCTCGGATGCCGTCTTA";
    for(size_t k:{11,23})for(bool reverse:{false,true}) {
        Graph g;g.k_=k;Options o;o.max_bp=900;o.max_coverage_penalty=100000;
        auto clean=dna(120)+observed_seed,sequence=clean+observed_tail+dna(150);
        auto spacer=observed_seed;
        if(reverse){sequence=reverse_complement(sequence);spacer=reverse_complement(spacer);clean=reverse_complement(clean);}
        g.add(sequence,12);Path p;
        for(size_t i=0;i+k<=spacer.size();++i)p.push_back(g.ids.at(spacer.substr(i,k)));
        auto seed=isolate(g,p,{},o).seeds.front();SearchDiagnostics diag;
        auto candidates=reconstruct(g,seed,{},o,nullptr,nullptr,nullptr,&diag);
        check(candidates.size()==1 && candidates[0].sequence==clean,"observed adapter branch rolled back to exact insert boundary on both strands");
        check(diag.adapter_stops==1 && diag.adapters[0].clean_flank_bp==0 && !diag.adapters[0].overlaps_seed,"boundary and contamination reason recorded");
        check((reverse?candidates[0].left:candidates[0].right).stop=="sequencing_adapter","adapter stop survives reconstruction");
        check((reverse?candidates[0].left:candidates[0].right).coverage.profile.samples==0,"removed adapter samples cannot improve ranking");
        o.adapters.clear();auto unsafe=reconstruct(g,seed,{},o);
        check(std::any_of(unsafe.begin(),unsafe.end(),[&](const auto& c){return c.sequence==sequence;}),"disabled guard reproduces original contaminated walk");
    }
    // Stable contaminating coverage must not beat a genuine continuation, even
    // when the valid path shares 12 adapter bases and diverges at the last one.
    Graph g;g.k_=23;Options o;o.max_bp=1000;o.max_coverage_penalty=100000;
    auto s=dna(34),left=dna(180),prefix=dna(150);
    auto real=left+s+prefix+adapter.substr(0,12)+"A"+dna(180);
    auto wrong=left+s+prefix+observed_tail+dna(180);
    g.add(real,12);g.add(wrong,12);Path p;
    for(auto& [ends,arc]:g.arcs)arc.count=12; // isolate adapter evidence at identical coverage
    for(size_t i=0;i+g.k()<=s.size();++i)p.push_back(g.ids.at(s.substr(i,g.k())));
    auto seed=isolate(g,p,{},o).seeds.front();SearchDiagnostics diag;
    auto candidates=reconstruct(g,seed,{},o,nullptr,nullptr,nullptr,&diag);
    check(std::any_of(candidates.begin(),candidates.end(),[&](const auto& c){return c.sequence==real;}),"genuine branch survives an almost-matching adapter prefix");
    check(std::none_of(candidates.begin(),candidates.end(),[&](const auto& c){return AdapterGuard(o.adapters).contains(c.sequence);}),"no reconstructed candidate crosses a complete adapter");
    check(diag.adapter_stops>0 && diag.adapters[0].clean_flank_bp==prefix.size(),"internal adapter boundary recorded in nucleotide coordinates");
    auto best=std::min_element(candidates.begin(),candidates.end(),[&](const auto& a,const auto& b){return better_prediction(a,s.size(),b,s.size(),o);});
    check(best->sequence==real,"clean full continuation wins over adapter-truncated alternative");
    // Adapter detection must include the seed/flank boundary, not just the
    // appended bases. It cannot silently relabel adapter bases as a protospacer.
    Graph overlap;overlap.k_=23;auto contaminated=dna(30)+adapter.substr(0,8);
    overlap.add(dna(100)+contaminated+observed_tail.substr(8)+dna(100));p.clear();
    for(size_t i=0;i+overlap.k()<=contaminated.size();++i)p.push_back(overlap.ids.at(contaminated.substr(i,overlap.k())));
    seed=isolate(overlap,p,{},o).seeds.front();diag={};
    candidates=reconstruct(overlap,seed,{},o,nullptr,nullptr,nullptr,&diag);
    check(candidates.size()==1 && candidates[0].sequence==contaminated && candidates[0].left.stop=="adapter_in_seed","seed-overlapping adapter rejected instead of retaining contaminated seed with clean opposite flank");
    check(diag.adapters.size()==1 && diag.adapters[0].overlaps_seed,"seed overlap diagnostic explicit");
    for(const auto& motif:o.adapters) {
        AdapterGuard guard(o.adapters);check(guard.contains("TT"+motif.sequence+"AA") && guard.contains(reverse_complement(motif.sequence)),"all default adapters checked on both strands");
    }
    check(!AdapterGuard(o.adapters).contains(adapter+"AAAAAA"),"isolated short core without adapter context is not enough to cut a biological path");
    check(AdapterGuard(o.adapters).contains(observed_tail),"observed adapter variant recognized with exact core and confirming bases");
    for(const auto& variant:{"AGATCGGAAGAGCCGTCGT", "AGATCGGAAGAGCACCACG", "AGATCGGAAGAGCGTCTGT", "AGATCGGAAGAGCACACT"}) {
        check(AdapterGuard(o.adapters).contains(variant),"observed adapter indel/short variant recognized");
        check(AdapterGuard(o.adapters).contains(reverse_complement(variant)),"adapter indel recognition is strand invariant");
    }
    Options exact;exact.adapters={{"AGATCGGAAGAGCACACGT"}};
    check(!AdapterGuard(exact.adapters).contains("AGATCGGAAGAGCACACT"),"user-specified adapters retain exact matching semantics");
    invalid([]{Options bad;bad.adapters={{"ACGT"}};bad.validate();});
    invalid([]{Options bad;bad.adapters={{"AGATCGGAAGANC"}};bad.validate();});
}
void test_graph_only_acceptance() {
    Fixture f;auto seed=f.full();auto all=reconstruct(f.g,seed,f.masked,f.o);
    auto it=std::find_if(all.begin(),all.end(),[&](const auto& c){return c.sequence==f.left+f.spacer+f.right;});
    check(it!=all.end(),"acceptance fixture has complete target");
    Candidate target=*it;target.external_left=target.left_length;target.external_right=target.right_length;
    auto original=refine_flanks(target,f.g.k(),f.o);
    check(original.size()==1 && original[0].sequence==target.sequence,"compatible supported flanks retained");
    auto profile=[&](double count){Coverage c;for(size_t i=0;i<10;++i)c.observe(count,f.o);return c;};
    target.left.coverage=profile(12);target.right.coverage=profile(100);
    auto split=refine_flanks(target,f.g.k(),f.o);
    check(split.size()==2 && split[0].flank_conflict && split[1].flank_conflict,"stable incompatible abundance contexts cannot be joined");
    check(split[0].sequence==f.left+f.spacer && split[1].sequence==f.spacer+f.right,"both one-sided contexts preserved as alternatives");
    for(const auto& c:split) {
        check(c.sequence.size()==c.path.size()+f.g.k()-1,"split flank keeps path-length accounting");
        check(c.sequence.substr(c.seed_begin,c.seed_end-c.seed_begin)==f.spacer,"split flank keeps seed coordinates");
        check((c.external_left==0)!=(c.external_right==0),"split candidate has exactly one external side");
    }
    Options legacy=f.o;legacy.joint_flank_check=false;
    check(refine_flanks(target,f.g.k(),legacy).size()==1,"joint compatibility check can be disabled explicitly");
    target.right.coverage=profile(14);
    check(refine_flanks(target,f.g.k(),f.o).size()==1,"small multiplicity differences are not a conflict");
    target.left.coverage={};target.left.coverage.observe(12,f.o);target.right.coverage=profile(100);
    check(refine_flanks(target,f.g.k(),f.o).size()==1,"too few observations cannot establish incompatible flank depth");
    Coverage dispersed;for(double count:{2.,20.,60.})dispersed.observe(count,f.o);
    check(joint_flank_loss(dispersed,profile(40),f.o)<f.o.max_coverage_penalty,"observed overdispersion is respected by joint check");
    // A short prefix with only seed-overlapping k-mers is not an independently
    // extended region. It remains unresolved if no other side can be retained.
    Candidate tiny;auto short_left=dna(f.g.k()),short_right=dna(f.g.k());
    tiny.sequence=short_left+f.spacer+short_right;tiny.path=f.g.add(tiny.sequence);
    tiny.seed_begin=short_left.size();tiny.seed_end=tiny.seed_begin+f.spacer.size();
    tiny.left_length=short_left.size();tiny.right_length=short_right.size();
    tiny.external_left=tiny.left_length;tiny.external_right=tiny.right_length;
    auto trimmed=refine_flanks(tiny,f.g.k(),f.o);
    check(trimmed.size()==1 && trimmed[0].sequence==f.spacer,"tails shorter than one complete external edge are removed");
    check(trimmed[0].left.context=="incomplete" && trimmed[0].right.context=="incomplete","short tails are unresolved, not a negative biological claim");
    check(trimmed[0].external_left==0 && trimmed[0].external_right==0 && trimmed[0].seed_begin==0,"short-tail trimming clears evidence and updates coordinates");
    legacy.min_external_flank=1;
    check(refine_flanks(tiny,f.g.k(),legacy)[0].sequence==tiny.sequence,"short-tail threshold is user adjustable");
    Candidate one=target;one.left_length=f.g.k();one.seed_begin=f.g.k();one.seed_end=one.seed_begin+f.spacer.size();
    one.sequence=short_left+f.spacer+f.right;one.path=f.g.add(one.sequence);one.left.coverage={};one.right.coverage=profile(12);
    auto kept=refine_flanks(one,f.g.k(),f.o);
    check(kept.size()==1 && kept[0].sequence==f.spacer+f.right,"weak side does not erase strong opposite flank");
    Candidate alternative=*it;
    alternative.sequence=dna(60)+f.left.substr(60)+f.spacer+f.right.substr(0,50)+dna(100);
    alternative.path=f.g.add(alternative.sequence);
    auto common=common_flanks(*it,0,alternative,0);
    check(common.seed_agrees && common.left>=40 && common.right>=50,"common sequence of competing paths is retained around the matched seed");
    Candidate consensus=*it;retain_flanks(consensus,common.left,common.right,f.g.k());
    check(f.g.spell(consensus.path)==consensus.sequence && consensus.sequence.substr(consensus.seed_begin,consensus.seed_end-consensus.seed_begin)==f.spacer,"ambiguity trimming preserves graph spelling and seed coordinates");
    alternative=*it;alternative.seed_begin+=5;alternative.seed_end-=3;
    common=common_flanks(*it,0,alternative,5);
    check(common.seed_agrees && common.left==it->seed_begin && common.right==it->right_length,"nested seed intervals compared in common query coordinates");
    alternative.sequence=dna(alternative.sequence.size());
    check(!common_flanks(*it,0,alternative,5).seed_agrees,"incompatible tied seed alignments cannot imply a shared locus");
}
int main() {
    try {
        test_extension();test_branches();test_coverage();test_reads_and_pam();test_join_and_repeated_seed();test_array_context();test_interval_compaction();test_stability_selection();test_adapter_branches();test_graph_only_acceptance();
        for(auto model:{"auto","poisson","negative-binomial"}){Options o;o.coverage_model=model;o.validate();}
        invalid([]{Options o;o.beam_width=0;o.validate();});
        invalid([]{Options o;o.max_coverage_penalty=std::nan("");o.validate();});
        invalid([]{Options o;o.coverage_model="bad";o.validate();});
        invalid([]{Options o;o.array_repeat_identity=1.1;o.validate();});
        std::cout<<checks<<" protospacer checks passed\n";
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
