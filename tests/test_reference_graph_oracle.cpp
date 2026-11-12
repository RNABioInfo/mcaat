#include "../bench/reference_graph_oracle.h"
#include "../bench/reference_root_graph.h"
#include "../bench/reference_unitig_oracle.h"
#include "../bench/reference_rolling_oracle.h"
#include "../bench/reference_segment_bounds.h"
#include "../bench/short_seed_placements.h"
#include "../bench/reachable_membership.h"
#include "../bench/adapter_filtered_graph.h"
#include <iostream>
#include <random>
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
};
// Independent string DP over every prefix of a single nonbranching path.
size_t optimal(const std::string& ref,const std::string& seq){std::vector<size_t> row(seq.size()+1),next(row.size());for(size_t j=0;j<row.size();++j)row[j]=j*9;for(size_t i=1;i<=ref.size();++i){next[0]=i*10;for(size_t j=1;j<row.size();++j)next[j]=std::min({row[j]+10,next[j-1]+9,row[j-1]+(ref[i-1]==seq[j-1]?0:10)});row.swap(next);}return *std::min_element(row.begin(),row.end());}
int main(){try{
    std::mt19937 rng(111129);auto dna=[&](size_t n){std::string s;while(n--)s+="ACGT"[rng()%4];return s;};
    for(bool back:{false,true}) {
        Graph base;auto sequence=dna(180);auto path=base.add(sequence);
        Node physical=back?path.at(130):path.at(10);auto word=base.label(physical);
        std::string prefix=back?word.substr(0,5):word.substr(18);if(back)std::reverse(prefix.begin(),prefix.end());
        audit::RootGraph<Graph> graph(base,back);Node start=graph.add_root(physical,prefix);
        audit::ReachableMembership<audit::RootGraph<Graph>> membership(graph,graph.native_membership_roots({start}),back);
        if(membership.find(word)==missing)throw std::runtime_error("cursor roots must not make their physical graph context appear disconnected");
        auto reference=back?sequence.substr(35,100):sequence.substr(28,100);if(back)std::reverse(reference.begin(),reference.end());
        auto found=audit::align(graph,start,reference,back,{100},100000);
        if(!found.found.count(100) || found.found.at(100).sequence!=reference)throw std::runtime_error("short-target cursor lost a real continuation");
    }
    size_t cases=0;
    for(bool back:{false,true}) {
        Graph graph;auto path=graph.add(dna(120));Node root=back?path.back():path.front();
        audit::RadiusMembership<Graph> radius(graph,{root},back,5);
        for(size_t i=0;i<path.size();++i) {
            auto node=back?path[path.size()-1-i]:path[i];
            if((radius.find(graph.label(node))!=missing)!=(i<=5))throw std::runtime_error("radius membership excludes a reachable node or includes an out-of-radius node");
        }
        audit::RadiusMembership<Graph> capped(graph,{root},back,5,1);
        if(capped.complete||capped.find(graph.label(path[10]))==missing)throw std::runtime_error("incomplete radius enumeration must retain global membership");
    }
    for(int variant=0;variant<24;++variant){auto seed=dna(23),ref=dna(100),seq=ref;
        if(variant<16)for(int m=0;m<variant;++m){size_t p=m*6;seq[p]="ACGT"[(std::string("ACGT").find(seq[p])+1)%4];}
        else if(variant<20)seq.insert(40,dna((variant-15)*3));else seq.erase(35,(variant-19)*3);
        Graph g;auto p=g.add(seed+seq);auto result=audit::align(g,p.front(),ref,false,{100},200000);
        bool expected=optimal(ref,seq)<=100;
        audit::SegmentBounds bounds;bounds.block=g.k();
        for(size_t p=0;p+g.k()<=ref.size();p+=g.k())bounds.costs.push_back(audit::segment_lower_bound(g,ref.substr(p,g.k())));
        auto pruned=audit::align(g,p.front(),ref,false,{100},200000,bounds);
        if(pruned.limited || bool(pruned.found.count(100))!=expected)throw std::runtime_error("segment bound changed oracle answer");
        std::set<size_t> lengths;bool any_expected=false;
        for(size_t length=40;length<=100;++length){lengths.insert(length);any_expected|=optimal(ref.substr(0,length),seq)<=length;}
        auto any=audit::align(g,p.front(),ref,false,lengths,200000,bounds,{},true);
        if(any.limited || !any.found.empty()!=any_expected)throw std::runtime_error("all-length search disagrees with independent DP");
        auto progress=audit::align(g,p.front(),ref,false,lengths,200000,bounds,{},true,true);
        if(progress.limited || !progress.found.empty()!=any_expected)throw std::runtime_error("progress ordering changed complete answer");
        audit::SegmentBounds stronger;stronger.block=g.k()+2;
        for(size_t pos=0;pos+stronger.block<=ref.size();pos+=stronger.block)
            stronger.costs.push_back(audit::two_edit_segment_lower_bound(g,ref.substr(pos,stronger.block)));
        stronger.prepare();
        auto lower=[&](size_t pos,size_t target){return std::max(bounds(pos,target),stronger(pos,target));};
        auto strong_result=audit::align(g,p.front(),ref,false,lengths,200000,lower,{},true);
        if(strong_result.limited || !strong_result.found.empty()!=any_expected)throw std::runtime_error("strong segment pruning changed complete answer");
        auto unitig_result=audit::align_unitigs(g,p.front(),ref,false,lengths,200000,lower,{},true);
        if(unitig_result.limited || !unitig_result.found.empty()!=any_expected)throw std::runtime_error("unitig alignment changed complete answer");
        auto rolling_result=audit::align_rolling(g,p.front(),ref,false,lengths,200000,lower,{},true);
        if(rolling_result.limited || !rolling_result.found.empty()!=any_expected)throw std::runtime_error("rolling alignment changed complete answer");
        audit::RadiusMembership<Graph> radius(g,{p.front()},false,112+g.k());
        audit::SegmentBounds radial;radial.block=g.k();
        for(size_t pos=0;pos+g.k()<=ref.size();pos+=g.k())radial.costs.push_back(audit::segment_lower_bound(radius,ref.substr(pos,g.k())));
        radial.prepare();
        auto radial_result=audit::align_rolling(g,p.front(),ref,false,lengths,200000,radial,{},true);
        if(radial_result.limited||!radial_result.found.empty()!=any_expected)throw std::runtime_error("radius bound changed independent DP result");
        auto rolling_limited=audit::align_rolling(g,p.front(),ref,false,{100},1);
        if(!rolling_limited.limited || rolling_limited.exhausted || rolling_limited.lowest_unexplored_cost>100)throw std::runtime_error("rolling cap produced a false absence certificate");
        if(result.found.count(100)!=expected || result.limited)throw std::runtime_error("oracle disagrees with independent string DP: "+std::to_string(variant));
        if(expected){auto w=result.found.at(100);if(w.matches*100<90*w.columns || seq.substr(0,w.sequence.size())!=w.sequence)throw std::runtime_error("bad witness");}
        auto backward_ref=ref;std::reverse(backward_ref.begin(),backward_ref.end());
        Graph bg;auto bp=bg.add(seq+seed);auto back=audit::align(bg,bp.back(),backward_ref,true,{100},200000);
        std::reverse(seq.begin(),seq.end());if(back.found.count(100)!=(optimal(backward_ref,seq)<=100))throw std::runtime_error("backward oracle mismatch");
        auto unitig_back=audit::align_unitigs(bg,bp.back(),backward_ref,true,{100},200000);
        if(unitig_back.limited || unitig_back.found.count(100)!=back.found.count(100))throw std::runtime_error("backward unitig mismatch");
        auto rolling_back=audit::align_rolling(bg,bp.back(),backward_ref,true,{100},200000);
        if(rolling_back.limited || rolling_back.found.count(100)!=back.found.count(100))throw std::runtime_error("backward rolling mismatch");
        std::reverse(seq.begin(),seq.end());
        auto bounded=audit::align(g,p.front(),ref,false,{24,100},120);
        if(!bounded.found.count(100) && bounded.lowest_unexplored_cost>100 && expected)
            throw std::runtime_error("false certificate of absent qualifying path");
        if(!bounded.found.count(24) && bounded.lowest_unexplored_cost>24 && optimal(ref.substr(0,24),seq)<=24)
            throw std::runtime_error("false short-prefix absence certificate");
        ++cases;
    }
    {
        Graph graph;auto seed=dna(23),adapter=Options{}.adapters.front().sequence;
        auto path=graph.add(seed+dna(12)+adapter+dna(50));
        std::unordered_map<Node,bool> cache;
        audit::AdapterFilteredGraph<Graph> filtered(graph,true,{path.front()},seed.size(),false,cache);
        bool stopped=false;for(size_t i=1;i<path.size();++i)if(!filtered.allowed(path[i])){stopped=true;break;}
        if(!stopped)throw std::runtime_error("adapter path not excluded");
        // A short seed beginning inside an adapter must not be rejected merely
        // because the un-emitted padding contains the rest of that adapter.
        auto label=adapter+"ACGTA";auto context=graph.add(label+dna(40));
        auto short_seed=label.substr(label.size()-17);
        audit::AdapterFilteredGraph<Graph> short_graph(graph,true,{context.front()},short_seed.size(),false,cache);
        for(size_t i=0;i<graph.k()-short_seed.size();++i)
            if(!short_graph.allowed(context[i]))throw std::runtime_error("short seed padding caused false exclusion");
    }
    // Every one-edit word needs a bound <= its cost, in either orientation.
    // A backward deletion can be a node prefix without a preceding node.
    for(size_t trial=0;trial<50;++trial) {
        auto ref=dna(23);
        for(size_t at=0;at<23;++at) {
            for(char b:std::string("ACGT")) {
                auto substituted=ref;substituted[at]=b;Graph sub;sub.add(substituted);
                if(audit::segment_lower_bound(sub,ref)>10)throw std::runtime_error("unsafe substitution bound");
                Graph ins;ins.add(ref.substr(0,at)+b+ref.substr(at));
                if(audit::segment_lower_bound(ins,ref)>9)throw std::runtime_error("unsafe insertion bound");
            }
            auto deleted=ref.substr(0,at)+ref.substr(at+1);
            Graph left,right;left.add("A"+deleted);right.add(deleted+"T");
            if(audit::segment_lower_bound(left,ref)>10 || audit::segment_lower_bound(right,ref)>10)
                throw std::runtime_error("unsafe deletion bound");
        }
    }
    // k+2 blocks cover two deletions without assuming an unobserved flanking
    // nucleotide. Check substitutions, insertions and both deletion ends.
    for(size_t trial=0;trial<8;++trial) {
        auto ref=dna(25);
        std::vector<std::string> variants{ref,ref.substr(2),ref.substr(0,23),"AC"+ref,
            ref.substr(0,9)+"AC"+ref.substr(9),ref.substr(0,10)+ref.substr(12)};
        auto changed=ref;changed[0]="ACGT"[(std::string("ACGT").find(changed[0])+1)%4];changed[24]="ACGT"[(std::string("ACGT").find(changed[24])+1)%4];variants.push_back(changed);
        for(auto seq:variants){Graph graph;graph.add(seq);if(audit::two_edit_segment_lower_bound(graph,ref)>20)throw std::runtime_error("unsafe two-edit block bound");}
    }
    Graph g;auto p=g.add(dna(160));auto small=audit::align(g,p.front(),dna(100),false,{100},2);
    if(!small.limited||small.exhausted)throw std::runtime_error("capped search mislabeled as absence");
    {
        Graph multi;auto wrong=multi.add(dna(130)),right=multi.add(dna(130));
        auto ref=multi.spell(right).substr(23,100);
        auto result=audit::align(multi,wrong.front(),ref,false,{100},200000,{},Path{right.front()});
        if(!result.found.count(100) || result.found.at(100).nodes.front()!=right.front())
            throw std::runtime_error("multi-start witness did not preserve its actual root");
        auto short_seed=multi.label(right.front()).substr(4);
        auto roots=audit::short_seed_endpoints(multi,short_seed,false);
        if(std::find(roots.begin(),roots.end(),right.front())==roots.end())throw std::runtime_error("short seed placement lost");
        for(auto node:roots)if(multi.label(node).substr(4)!=short_seed)throw std::runtime_error("wrong short seed placement");
        audit::ReachableMembership<Graph> membership(multi,Path{right.front()},false);
        if(membership.find(multi.label(wrong.back()))!=missing)throw std::runtime_error("disconnected component accepted as reachable");
        for(auto node:right)if(membership.find(multi.label(node))!=node)throw std::runtime_error("reachable member rejected");
        audit::ReachableMembership<Graph> backward(multi,Path{right.back()},true);
        for(auto node:right)if(backward.find(multi.label(node))!=node)throw std::runtime_error("backward reachable member rejected");
        audit::ReachableMembership<Graph> limited(multi,Path{right.front()},false);limited.work_limit=2;
        if(limited.find(multi.label(wrong.back()))==missing || !limited.unknown)throw std::runtime_error("bounded reachability produced a false negative");
    }
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
        auto compact=audit::align_unitigs(graph,path.front(),ref,false,{10},50000);
        if(compact.limited || bool(compact.found.count(10))!=expected)throw std::runtime_error("unitig oracle disagrees with exhaustive walks");
        auto rolling=audit::align_rolling(graph,path.front(),ref,false,{10},50000);
        if(rolling.limited || bool(rolling.found.count(10))!=expected)throw std::runtime_error("rolling oracle disagrees with exhaustive walks");
        std::vector<FlankContext> models{{"A",native.substr(4,3),""}};
        bool permitted=false;
        std::function<void(Node,std::string)> filtered_visit=[&](Node node,std::string sequence) {
            auto context=context_state(sequence,models,1.0);
            if(context==ContextState::array)return;
            if(context==ContextState::external && optimal(ref,sequence)<=ref.size()){permitted=true;return;}
            if(sequence.size()>=11 || permitted)return;
            for(auto arc:graph.neighbors(node,false))filtered_visit(arc.to,sequence+graph.extension_base(arc,false));
        };
        filtered_visit(path.front(),"");
        audit::PrefixContext prefix(models,1.0);
        auto filtered=audit::align(graph,path.front(),ref,false,{10},50000,{}, {},false,false,&prefix);
        if(filtered.limited || bool(filtered.found.count(10))!=permitted)
            throw std::runtime_error("context-aware oracle disagrees with exhaustive valid graph walks");
        audit::PrefixContext compact_prefix(models,1.0);
        auto compact_filtered=audit::align_unitigs(graph,path.front(),ref,false,{10},50000,{}, {},false,false,&compact_prefix);
        if(compact_filtered.limited || bool(compact_filtered.found.count(10))!=permitted)throw std::runtime_error("unitig source context disagrees with exhaustive walks");
        audit::PrefixContext rolling_prefix(models,1.0);
        auto rolling_filtered=audit::align_rolling(graph,path.front(),ref,false,{10},50000,{}, {},false,&rolling_prefix);
        if(rolling_filtered.limited || bool(rolling_filtered.found.count(10))!=permitted)throw std::runtime_error("rolling source context disagrees with exhaustive walks");
    }
    {
        Graph graph;graph.k_=3;graph.add("AACGTAACGTAACGT");graph.add("AACGGGAAACGT");
        std::vector<FlankContext> models{{"","","AACGT"}};
        for(bool back:{false,true})for(auto ref:std::vector<std::string>{"GTAACGTAAC","GGGAAACGTA","CAATGCAATG"}) {
            auto root=graph.find("AAC");SourceCycleGuard guard(graph.k(),&models);
            audit::CycleContext automaton(graph.k(),models);bool expected=false;
            std::function<void(Path,size_t,size_t,std::string)> visit=[&](Path nodes,size_t run,size_t state,std::string sequence) {
                if(optimal(ref,sequence)<=ref.size())expected=true;
                if(sequence.size()>=11)return;
                Side side;side.nodes=nodes;
                for(auto arc:graph.neighbors(nodes.back(),back)) {
                    auto next_run=guard.contains(graph,nodes.back(),arc,back)?run+1:0;
                    bool rejected=!repeat_stop(Path{root},side,arc.to,next_run,64).empty();
                    auto next_state=automaton.advance(graph,state,nodes.back(),arc,back);
                    if(rejected!=(next_state==audit::CycleContext::reject))throw std::runtime_error("source cycle automaton changed the production rule");
                    if(rejected)continue;
                    auto next=nodes;next.push_back(arc.to);visit(next,next_run,next_state,sequence+graph.extension_base(arc,back));
                }
            };
            visit({root},0,0,"");audit::CycleContext alignment_context(graph.k(),models);
            auto result=audit::align(graph,root,ref,back,{10},200000,{}, {},false,false,nullptr,&alignment_context);
            if(result.limited || bool(result.found.count(10))!=expected)throw std::runtime_error("cycle-aware alignment disagrees with exhaustive walks");
            audit::CycleContext compact_cycles(graph.k(),models);
            auto compact=audit::align_unitigs(graph,root,ref,back,{10},200000,{}, {},false,false,nullptr,&compact_cycles);
            if(compact.limited || bool(compact.found.count(10))!=expected)throw std::runtime_error("unitig cycle-aware alignment disagrees with exhaustive walks");
            audit::CycleContext rolling_cycles(graph.k(),models);
            auto rolling=audit::align_rolling(graph,root,ref,back,{10},200000,{}, {},false,nullptr,&rolling_cycles);
            if(rolling.limited || bool(rolling.found.count(10))!=expected)throw std::runtime_error("rolling cycle-aware alignment disagrees with exhaustive walks");
        }
    }
    std::cout<<cases<<" forward/backward DP cases and 12 exhaustive branching graphs passed; capped search stays unresolved\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
