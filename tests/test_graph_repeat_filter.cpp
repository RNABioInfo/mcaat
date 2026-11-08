#include "crispr_array/graph_repeat_filter.h"
#include <iostream>
#include <random>
#include <sstream>
using namespace mcaat::arrays;
size_t checks=0;
void check(bool ok,const std::string& message) {
    ++checks;if(!ok) throw std::runtime_error(message);
}
struct Graph {
    size_t k;
    std::map<std::string,Node> ids;
    std::map<Node,std::string> labels;
    std::set<std::pair<Node,Node>> links;
    Node node(const std::string& s) {
        if (!ids.count(s)) {Node id=ids.size()+1;ids[s]=id;labels[id]=s;}
        return ids.at(s);
    }
    Path add_one(const std::string& period) {
        Path path;auto padded=period+period.substr(0,k-1);
        for(size_t i=0;i<period.size();++i) path.push_back(node(padded.substr(i,k)));
        for(size_t i=0;i<path.size();++i) links.emplace(path[i],path[(i+1)%path.size()]);
        return path;
    }
    Path add(const std::string& period) {
        auto p=add_one(period);add_one(reverse_complement(period));return p;
    }
    std::string label(Node n) {return labels.at(n);}
    bool adjacent(Node a,Node b) {return links.count({a,b});}
    Node reverse_node(Node n) {
        auto it=ids.find(reverse_complement(label(n)));
        return it==ids.end()?missing_node:it->second;
    }
};
std::mt19937 rng(107);
std::string dna(size_t n) {std::string r;while(n--) r+="ACGT"[rng()%4];return r;}
std::string mutate(std::string r,size_t pos) {r[pos]=r[pos]=='A'?'C':'A';return r;}
struct Fixture {
    Graph g;
    std::string repeat;
    std::vector<Path> cycles;
    std::vector<std::string> spacers,repeats;
    Options o;
    Fixture(size_t k=9,size_t r=28):g{k},repeat(dna(r)) {
        for (;;) {
            bool valid=true;
            for (auto rep:{repeat,mutate(repeat,0),mutate(repeat,std::min<size_t>(12,r-1)),mutate(repeat,std::min<size_t>(14,r-1)),mutate(repeat,std::min<size_t>(26,r-1))}) {
                std::set<std::string> words;
                for (size_t i=0;i+k<=rep.size();++i) words.insert(rep.substr(i,k));
                valid &= words.size()==r-k+1;
            }
            if (valid) break;
            repeat=dna(r);
        }
        repeat=std::min(repeat,reverse_complement(repeat));
        o.k=k;o.repeat_min=r;o.repeat_max=r;o.spacer_min=20;o.spacer_max=50;
    }
    void append(const std::string& rep,size_t length=30) {
        for(size_t attempt=0;attempt<10000;++attempt) {
            auto s=dna(length);s.front()="ACGT"[(cycles.size()+attempt)%4];s.back()="TGCA"[(cycles.size()+attempt)%4];
            auto period=rep+s,full=period+period.substr(0,g.k-1);
            std::set<std::string> words;
            for(size_t i=0;i<period.size();++i) words.insert(full.substr(i,g.k));
            if(words.size()!=period.size()) continue;
            cycles.push_back(g.add(period));spacers.push_back(s);repeats.push_back(rep);return;
        }
        throw std::runtime_error("fixture generation failed: repeat="+rep+" k="+std::to_string(g.k)+" spacer="+std::to_string(length));
    }
    std::vector<System> run(std::map<Node,std::vector<Path>> groups={}) {
        RepeatFilter<Graph> f(g,o);
        return f.run(groups.empty()?std::map<Node,std::vector<Path>>{{1,cycles}}:groups);
    }
};
void test_threshold_and_variants() {
    Fixture f;
    for(int i=0;i<20;++i) f.append(i<17?f.repeat:mutate(f.repeat,14));
    auto systems=f.run();check(systems.size()==1,"17/20 default accepted");
    auto& s=systems.front();check(s.support_count==17,"support numerator is distinct exact full paths");
    check(s.candidate_cycles.size()==20 && s.observations.size()==20,"anchored minority variants retained");
    std::set<std::pair<std::string,std::string>> actual,expected;
    for(auto& row:s.observations) {
        actual.emplace(row.repeat,row.spacer);
        auto parts=partition(row.cycle,row.repeat_nodes,f.o.k);
        check(parts.repeat.size()==28-9+1 && parts.spacer.size()==30-9+1,"pure k-mer parts");
        check(parts.repeat_to_spacer.size()==8 && parts.spacer_to_repeat.size()==8,"boundary k-mers preserved");
        std::set<Node> all;
        for(auto p:{parts.repeat,parts.spacer,parts.repeat_to_spacer,parts.spacer_to_repeat}) all.insert(p.begin(),p.end());
        check(all.size()==row.cycle.size(),"partition covers whole cycle");
    }
    for(size_t i=0;i<20;++i) expected.emplace(f.repeats[i],f.spacers[i]);
    check(actual==expected,"observed sequences preserved exactly, including variants");
    f.o.repeat_support=.851;check(f.run().empty(),"threshold is inclusive but not rounded down");
    f.o.repeat_support=.85;
    auto copies=f.cycles;
    for(int i=0;i<50;++i) copies.push_back(rotate(f.cycles[0],i%f.cycles[0].size()));
    RepeatFilter<Graph> filter(f.g,f.o);systems=filter.run({{1,copies}});
    check(systems.size()==1 && systems[0].support_count==17 && filter.stats.duplicate_cycles==50,"duplicates cannot inflate support");
    Fixture sixty;
    for(int i=0;i<20;++i) sixty.append(i<12?sixty.repeat:mutate(sixty.repeat,14));
    check(sixty.run().empty(),"12/20 rejected at .85");
    sixty.o.repeat_support=.6;systems=sixty.run();
    check(systems.size()==1 && systems[0].support_count==12,"12/20 accepted at .60");
}
void test_joint_support() {
    Fixture f(9,40);
    for(int i=0;i<20;++i) f.append(i<3?mutate(f.repeat,12):(i<6?mutate(f.repeat,26):f.repeat));
    check(f.run().empty(),"separate 85% node majorities cannot create an 85% whole repeat");
    f.o.repeat_support=.7;auto systems=f.run();
    check(systems.size()==1 && systems[0].support_count==14,"observed whole path accepted at actual support");
}
void test_deduplication() {
    Fixture f;
    for(int i=0;i<8;++i) f.append(f.repeat);
    auto a=std::vector<Path>(f.cycles.begin(),f.cycles.begin()+4);
    auto b=std::vector<Path>(f.cycles.begin()+2,f.cycles.begin()+6);
    auto c=std::vector<Path>(f.cycles.begin()+4,f.cycles.end());
    RepeatFilter<Graph> filter(f.g,f.o);
    auto systems=filter.run({{1,a},{2,b}});
    check(systems.size()==1 && systems[0].observations.size()==6 && filter.stats.merged_systems==1,"overlapping systems merge and add spacers");
    systems=filter.run({{1,a},{2,f.cycles}});
    check(systems.size()==1 && systems[0].observations.size()==8 && filter.stats.duplicate_systems==1,"subset detection removed");
    check(systems[0].sources.size()==2,"source starts retained");
    systems=filter.run({{1,a},{2,c}});
    check(systems.size()==2,"shared repeat with disjoint cycles not merged");
    auto rotated=a;
    for(auto& p:rotated) p=rotate(p,17);
    systems=filter.run({{1,a},{2,rotated}});
    check(systems.size()==1,"different start positions deduplicate");
    auto rc=a;
    for(auto& p:rc) {std::reverse(p.begin(),p.end());for(auto& n:p) n=f.g.reverse_node(n);}
    systems=filter.run({{1,a},{2,rc}});
    check(systems.size()==1,"reverse-complement detections deduplicate");
    std::vector<Path> different_repeat;
    auto other=dna(f.repeat.size());
    for(size_t i=0;i<4;++i) different_repeat.push_back(f.g.add(other+f.spacers[i]));
    systems=filter.run({{1,a},{2,different_repeat}});
    check(systems.size()==2,"shared spacers with distinct repeat paths remain separate");
}
void test_shifted_strands() {
    Graph graph{23};
    const std::string repeat="GTTTCAATCCACGCGCCCGCACGGGGCGCGAC";
    std::vector<Path> cycles;
    for(size_t i=0;i<8;++i) {
        auto spacer=dna(35);spacer.front()="ACGT"[i%4];
        spacer[33]="ACGT"[i%4];spacer[34]=i<4?'A':"CGTC"[i%4];
        cycles.push_back(graph.add(repeat+spacer));
    }
    Options o;o.k=23;o.repeat_min=repeat.size();o.repeat_max=repeat.size()+1;
    RepeatFilter<Graph> filter(graph,o);
    std::vector<Path> first(cycles.begin(),cycles.begin()+4);
    auto individual=filter.run({{1,first}});
    check(individual.size()==1 && filter.spell(individual[0].repeat_path)=="A"+repeat,"left extension changes canonical strand");
    auto all=filter.run({{1,first},{2,cycles}});
    check(all.size()==1 && all[0].observations.size()==8,"shifted repeat boundaries and opposite canonical strands still deduplicate");
}
void test_incremental_deduplication() {
    Fixture f(23,30);
    for (int i=0;i<16;++i) f.append(f.repeat);
    std::map<Node,std::vector<Path>> groups;
    const size_t copies=256;
    for (size_t i=0;i<copies;++i) groups[i]=f.cycles;
    RepeatFilter<Graph> filter(f.g,f.o);
    auto systems=filter.run(groups);
    check(systems.size()==1 && systems[0].sources.size()==copies,"all duplicate starts retained");
    check(filter.stats.duplicate_systems==copies-1,"all duplicate systems removed");
    check(filter.stats.dedup_index_builds==copies,"source-only removals never rebuild system index");
    check(filter.stats.dedup_pair_checks==copies-1,"dense duplicate groups require only linear pair checks");
    // Each superset replaces the retained system's geometry and index. Old
    // entries must disappear, and every source remains attached to the union.
    groups.clear();
    for (size_t i=4;i<=16;i+=4) groups[i]={f.cycles.begin(),f.cycles.begin()+i};
    systems=filter.run(groups);
    check(systems.size()==1 && systems[0].observations.size()==16,"successive supersets replace indexed geometry");
    check(systems[0].sources==Path({4,8,12,16}),"superset replacements preserve sorted provenance");
    check(filter.stats.dedup_index_builds==7,"only changed superset geometry gets reindexed");
    // Disjoint groups become connected through a later candidate.
    groups={{1,{f.cycles.begin(),f.cycles.begin()+4}},
            {2,{f.cycles.begin()+4,f.cycles.begin()+8}},
            {3,{f.cycles.begin()+2,f.cycles.begin()+6}}};
    systems=filter.run(groups);
    check(systems.size()==1 && systems[0].observations.size()==8,"bridge merge discovers newly adjacent group");
    check(systems[0].sources==Path({1,2,3}) && filter.stats.merged_systems==2,"bridge merge preserves sources and counts");
}
void test_boundaries_and_rejections() {
    Fixture f;
    for(int i=0;i<20;++i) f.append(i<17?f.repeat:mutate(f.repeat,0));
    auto systems=f.run();
    check(systems.size()==1 && systems[0].candidate_cycles.size()==20 && systems[0].observations.size()==17,"missing boundary anchors remain unresolved");
    Fixture one;
    one.append(one.repeat);check(one.run().empty(),"single cycle rejected");
    auto invalid=one.cycles[0];invalid.push_back(invalid[2]);
    RepeatFilter<Graph> filter(one.g,one.o);
    systems=filter.run({{1,{invalid,one.cycles[0]}}});
    check(systems.empty() && filter.stats.invalid_cycles==1,"non-simple or disconnected cycle rejected");
    systems=filter.run({{1,{invalid,invalid,one.cycles[0]}}});
    check(systems.empty() && filter.stats.invalid_cycles==2 && filter.stats.duplicate_cycles==0,
          "cached invalid copies still count as invalid input cycles");
    Fixture cached;
    for (int i=0;i<4;++i) cached.append(cached.repeat);
    RepeatFilter<Graph> reused(cached.g,cached.o);
    check(reused.run({{1,cached.cycles}}).size()==1,"validation cache populated");
    cached.g.links.erase({cached.cycles[0][35],cached.cycles[0][36]});
    reused.run({{1,cached.cycles}});
    check(reused.stats.invalid_cycles==1,"validation cache cleared between runs after graph changes");
    Fixture k_repeat(9,9);
    for(int i=0;i<4;++i) k_repeat.append(k_repeat.repeat);
    systems=k_repeat.run();check(systems.size()==1 && systems[0].repeat_path.size()==1,"a repeat exactly k long is a shared vertex");
    Fixture short_spacer(9,12);short_spacer.o.spacer_min=6;
    for(int i=0;i<4;++i) short_spacer.append(short_spacer.repeat,6);
    systems=short_spacer.run();check(systems.size()==1,"spacer shorter than k supported without fabricated pure spacer nodes");
    for(const auto& row:systems[0].observations) check(partition(row.cycle,row.repeat_nodes,9).spacer.empty(),"short spacer has boundary nodes only");
}
System assigned(Fixture& f,const std::vector<size_t>& indexes) {
    System s;s.repeat_path=Path(f.cycles[0].begin(),f.cycles[0].begin()+f.repeat.size()-f.g.k+1);
    for(auto i:indexes) {
        s.candidate_cycles.push_back(cycle_key(f.cycles.at(i)));
        s.observations.push_back({f.cycles[i],f.repeat.size()-f.g.k+1,f.repeats[i],f.spacers[i],true});
    }
    s.support_count=s.candidate_count=indexes.size();return s;
}
void test_curation_ownership() {
    Fixture f(23,28);f.o.quality.strict=true;for(int i=0;i<10;++i)f.append(f.repeat);
    auto a=assigned(f,{0,1,2,3,4,5}),b=assigned(f,{0,1,6,7});
    RepeatFilter<Graph> filter(f.g,f.o);
    auto systems=filter.curate({a,b});
    check(systems.size()==2 && systems[0].observations.size()==6 && systems[1].observations.size()==2,"independent residual evidence survives shared-cycle removal");
    check(systems[1].support_count==2 && systems[1].candidate_count==2,"support recomputed on residual evidence");
    check(filter.decisions[1].removed_cycles==2 && filter.decisions[1].competing_candidates.count(1),"competing candidate provenance recorded");
    systems=filter.curate({a,assigned(f,{0,1,6})});
    check(systems.size()==1 && filter.decisions[1].status=="shared_cycle_subset","one residual spacer cannot remain an array");
    auto ambiguous=b;ambiguous.candidate_count=5;
    systems=filter.curate({a,ambiguous});
    check(systems.size()==1 && filter.decisions[1].status=="repeat_support_after_overlap","unresolved cycles stay in support denominator");
    auto bad=a;bad.candidate_count=500;
    systems=filter.curate({bad,b});
    check(systems.size()==1 && systems[0].observations.size()==4 && filter.decisions[0].status=="cycle_limit","rejected candidates never claim cycles");
    auto reverse=a;reverse.repeat_path.clear();reverse.observations.clear();reverse.candidate_cycles.clear();
    for(const auto& r:a.observations) {
        auto repeat=reverse_complement(r.repeat),spacer=reverse_complement(r.spacer);auto c=f.g.add(repeat+spacer);
        reverse.observations.push_back({c,r.repeat_nodes,repeat,spacer,true});reverse.candidate_cycles.push_back(rotate(c,13));
        if(reverse.repeat_path.empty())reverse.repeat_path=Path(c.begin(),c.begin()+r.repeat_nodes);
    }
    systems=filter.curate({a,reverse});
    check(systems.size()==1 && systems[0].observations.size()==6,"rotated reverse-strand circles have one accepted owner");
    Fixture pal(9,28);auto half=dna(14);pal.repeat=half+reverse_complement(half);
    for(int i=0;i<3;++i)pal.append(pal.repeat);
    auto both=assigned(pal,{0,1,2});
    for(size_t i=0;i<3;++i) {
        auto spacer=reverse_complement(pal.spacers[i]);auto cycle=pal.g.add(pal.repeat+spacer);
        both.candidate_cycles.push_back(cycle_key(cycle));
        both.observations.push_back({cycle,20,pal.repeat,spacer,true});
    }
    both.support_count=both.candidate_count=6;RepeatFilter<Graph> pal_filter(pal.g,pal.o);
    auto pal_systems=pal_filter.curate({both});
    check(pal_systems.size()==1 && pal_systems[0].observations.size()==3,"opposite strand copies inside one candidate cannot double-count cycles");
    check(pal_systems[0].candidate_count==3 && pal_systems[0].support_count==3 && pal_filter.decisions[0].duplicate_strand_cycles==3,"within-candidate strand duplicates adjust support denominator and remain traceable");
    auto conservative_options=f.o;conservative_options.quality.strict=false;
    RepeatFilter<Graph> conservative(f.g,conservative_options);
    auto preserved=conservative.curate({a,b});
    check(preserved.size()==2 && preserved[0].observations.size()==6 && preserved[1].observations.size()==4,"shared evidence must not be stolen from another candidate in conservative mode");
    check(conservative.decisions[0].status=="accepted_with_shared_evidence" && conservative.decisions[1].status=="accepted_with_shared_evidence","shared evidence must be flagged symmetrically");
    check(preserved[0].quality.flags.back()=="shared_cycle_evidence" && preserved[1].quality.flags.back()=="shared_cycle_evidence","ambiguous overlap is flagged on both candidates");
    conservative_options.cycle_min=59;conservative_options.cycle_max=100;
    RepeatFilter<Graph> outside(f.g,conservative_options);
    check(outside.curate({a}).empty() && outside.decisions[0].status=="cycle_length_range","saved candidates outside requested cycle bounds are rejected explicitly");
    // Shift an otherwise identical circle so the proposed repeat is inside the
    // original spacer. It is a competing boundary interpretation, not new data.
    auto wrong=assigned(f,{0,1});wrong.repeat_path.clear();
    for(auto& row:wrong.observations) {
        auto period=row.repeat+row.spacer;size_t offset=row.repeat.size();period=period.substr(offset)+period.substr(0,offset);
        row.cycle=rotate(row.cycle,offset);row.repeat=period.substr(0,28);row.spacer=period.substr(28);
        if(wrong.repeat_path.empty())wrong.repeat_path=Path(row.cycle.begin(),row.cycle.begin()+row.repeat_nodes);
    }
    systems=filter.curate({a,wrong});
    check(systems.size()==1 && filter.decisions[1].retained_paths.empty(),"repeat/spacer reinterpretation cannot duplicate accepted evidence");
    check(filter.decisions[1].original.observations.size()==2,"alternative boundary observations preserved for review");
    for(const auto& row:systems[0].observations)check(row.repeat==f.repeat,"winning repeat boundary unchanged");
}

std::set<std::string> circles(const std::vector<System>& systems) {
    std::set<std::string> result;
    for(const auto& s:systems)for(const auto& row:s.observations)
        result.insert(circular_sequence_key(row.repeat+row.spacer));
    return result;
}
System shift_boundary(System s,size_t offset,size_t length,size_t k) {
    s.repeat_path.clear();
    for(auto& row:s.observations) {
        auto sequence=row.repeat+row.spacer;
        sequence=sequence.substr(offset)+sequence.substr(0,offset);
        row.cycle=rotate(row.cycle,offset);row.repeat=sequence.substr(0,length);row.spacer=sequence.substr(length);
        row.repeat_nodes=length-k+1;
        if(s.repeat_path.empty())s.repeat_path=Path(row.cycle.begin(),row.cycle.begin()+row.repeat_nodes);
    }
    return s;
}
void test_sequence_subset_removal() {
    Fixture f(23,28);for(int i=0;i<10;++i)f.append(f.repeat);
    auto options=f.o;options.repeat_min=23;
    RepeatFilter<Graph> filter(f.g,options);
    auto parent=assigned(f,{0,1,2,3,4,5});parent.sources={10};
    auto child=shift_boundary(assigned(f,{0,1,2}),2,26,23);child.sources={20};
    const auto before=circles({parent,child});
    auto result=filter.curate({child,parent});
    check(result.size()==1 && result.front().observations.size()==6,"contained circles removed despite different boundary and spacers");
    check(circles(result)==before,"subset removal preserves every distinct accepted circular sequence");
    check(filter.decisions[0].status=="redundant_circular_subset" && filter.decisions[0].competing_candidates.count(2),"redundant report links to a retained parent");
    check(filter.decisions[0].retained_paths.empty() && filter.decisions[0].original.observations.size()==3,"removed boundary alternatives remain traceable");
    check(result[0].sources==std::vector<Node>({10,20}),"source provenance merged without increasing support");
    check(filter.stats.redundant_candidates==1 && filter.stats.redundant_observations==3,"redundancy counts are explicit");
    result=filter.curate(result);
    check(result.size()==1 && filter.stats.redundant_candidates==0 && circles(result)==before,"sequence deduplication is idempotent");

    auto equal=shift_boundary(parent,2,26,23);
    for(auto input:{std::vector<System>{equal,parent},std::vector<System>{parent,equal}}) {
        result=filter.curate(input);
        check(result.size()==1 && filter.spell(result[0].repeat_path)==f.repeat,"equal evidence picks the same existing eligible boundary independent of input order");
    }
    auto reverse=child;reverse.repeat_path.clear();reverse.candidate_cycles.clear();
    for(auto& row:reverse.observations) {
        row.repeat=reverse_complement(row.repeat);row.spacer=reverse_complement(row.spacer);
        row.cycle=f.g.add(row.repeat+row.spacer);reverse.candidate_cycles.push_back(cycle_key(row.cycle));
        if(reverse.repeat_path.empty())reverse.repeat_path=Path(row.cycle.begin(),row.cycle.begin()+row.repeat_nodes);
    }
    result=filter.curate({reverse,parent});
    check(result.size()==1 && circles(result)==before,"subset matching handles circular rotation and reverse complement together");

    result=filter.curate({assigned(f,{0,1}),assigned(f,{1,2}),assigned(f,{0,2})});
    check(result.size()==3,"coverage by a union of different candidates never causes deletion");
    result=filter.curate({parent,assigned(f,{0,1,6})});
    check(result.size()==2 && circles(result)==circles({parent,assigned(f,{0,1,6})}),"one unique circular sequence prevents subset removal");
    result=filter.curate({assigned(f,{0,1}),assigned(f,{6,7})});
    check(result.size()==2,"identical repeats with distinct circles are not duplicate systems");

    auto unresolved=child;unresolved.candidate_count=4;unresolved.candidate_cycles.push_back(cycle_key(f.cycles[9]));
    options.repeat_support=.6;RepeatFilter<Graph> allow_unresolved(f.g,options);
    result=allow_unresolved.curate({unresolved,parent});
    check(result.size()==1 && allow_unresolved.decisions[0].original.candidate_cycles.size()==4,
          "unassigned source cycles survive in the removed candidate's review evidence");

    std::vector<System> duplicates(256,parent);
    result=filter.curate(std::move(duplicates));
    check(result.size()==1 && filter.stats.sequence_subset_checks==255,"identical reports require linear subset checks, not a pairwise clique");
}
void test_curation_length_bounds() {
    Fixture f(23,28);for(size_t n:{29,30,30,31})f.append(f.repeat,n);
    auto options=f.o;options.repeat_min=options.repeat_max=28;
    options.spacer_min=options.spacer_max=30;options.repeat_support=.5;
    auto input=assigned(f,{0,1,2,3});
    RepeatFilter<Graph> filter(f.g,options);
    auto result=filter.curate({input});
    check(result.size()==1 && result[0].observations.size()==2,"curation enforces inclusive spacer bounds, not just repeat inference");
    check(result[0].candidate_count==4 && result[0].support_count==2,"invalid length assignments do not erase the original support denominator");
    check(filter.stats.invalid_length_observations==2 && filter.decisions[0].invalid_length_observations==2,"length exclusions are counted");
    for(const auto& row:result[0].observations)check(row.repeat.size()==28 && row.spacer.size()==30,"boundary lengths retained exactly");
    options.repeat_min=29;options.repeat_max=30;RepeatFilter<Graph> too_short(f.g,options);
    check(too_short.curate({input}).empty() && too_short.decisions[0].status=="repeat_length_range","below-minimum repeat rejected in curation");
    options.repeat_min=23;options.repeat_max=27;RepeatFilter<Graph> too_long(f.g,options);
    check(too_long.curate({input}).empty() && too_long.decisions[0].status=="repeat_length_range","above-maximum repeat rejected in curation");
    options.repeat_max=28;options.spacer_min=options.spacer_max=32;RepeatFilter<Graph> no_spacers(f.g,options);
    check(no_spacers.curate({input}).empty() && no_spacers.decisions[0].status=="repeat_spacer_length_range","candidate with no length-valid spacers rejected explicitly");
}
void test_explicit_string_exclusions() {
    Fixture f(23,28);for(int i=0;i<6;++i)f.append(f.repeat);
    RepeatFilter<Graph> filter(f.g,f.o);
    auto parent=assigned(f,{0,1,2,3,4,5}),child=assigned(f,{0,1,2});
    auto result=filter.curate({parent,child},{},{{1,"short_motif_across_unit"}});
    check(result.size()==1 && result[0].candidate_id==2 && result[0].observations.size()==3,
          "explicit string rejection occurs before containment and cannot suppress an unflagged child");
    check(filter.decisions[0].status=="string_filter:short_motif_across_unit" &&
          filter.decisions[0].original.observations.size()==6 && filter.decisions[0].retained_paths.empty(),
          "string-rejected observations and exact reason remain in the review decisions");
    check(filter.stats.quality_rejected==1 && filter.stats.redundant_candidates==0,
          "string rejection and redundancy are counted separately");
    for(size_t id:{0,3}) {
        bool failed=false;
        try {filter.curate({parent,child},{},{{id,"short_motif_across_unit"}});}
        catch(const std::invalid_argument&) {failed=true;}
        check(failed,"unknown explicit exclusion ID rejected");
    }
}

int main() {
    try {
        for(auto s:{"0","-0.1","1.01","nan","inf","85%",".6junk",""}) {
            bool threw=false;try {parse_repeat_support(s);}catch(const std::invalid_argument&){threw=true;}
            check(threw,std::string("invalid support rejected: ")+s);
        }
        check(parse_repeat_support("0.85")==.85 && parse_repeat_support(".60")==.6 && parse_repeat_support("1")==1,"valid fractions");
        test_threshold_and_variants();test_joint_support();test_deduplication();test_shifted_strands();
        test_incremental_deduplication();test_boundaries_and_rejections();test_curation_ownership();
        test_sequence_subset_removal();test_curation_length_bounds();
        test_explicit_string_exclusions();
        std::cout<<checks<<" graph repeat-filter checks passed\n";
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
