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
    Fixture k_repeat(9,9);
    for(int i=0;i<4;++i) k_repeat.append(k_repeat.repeat);
    systems=k_repeat.run();check(systems.size()==1 && systems[0].repeat_path.size()==1,"a repeat exactly k long is a shared vertex");
    Fixture short_spacer(9,12);short_spacer.o.spacer_min=6;
    for(int i=0;i<4;++i) short_spacer.append(short_spacer.repeat,6);
    systems=short_spacer.run();check(systems.size()==1,"spacer shorter than k supported without fabricated pure spacer nodes");
    for(const auto& row:systems[0].observations) check(partition(row.cycle,row.repeat_nodes,9).spacer.empty(),"short spacer has boundary nodes only");
}
int main() {
    try {
        for(auto s:{"0","-0.1","1.01","nan","inf","85%",".6junk",""}) {
            bool threw=false;try {parse_repeat_support(s);}catch(const std::invalid_argument&){threw=true;}
            check(threw,std::string("invalid support rejected: ")+s);
        }
        check(parse_repeat_support("0.85")==.85 && parse_repeat_support(".60")==.6 && parse_repeat_support("1")==1,"valid fractions");
        test_threshold_and_variants();test_joint_support();test_deduplication();test_shifted_strands();
        test_incremental_deduplication();test_boundaries_and_rejections();
        std::cout<<checks<<" graph repeat-filter checks passed\n";
    }catch(const std::exception& e){std::cerr<<"FAILED: "<<e.what()<<'\n';return 1;}
}
