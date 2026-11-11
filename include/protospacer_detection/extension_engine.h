#ifndef MCAAT_PROTO_EXTENSION_ENGINE_H
#define MCAAT_PROTO_EXTENSION_ENGINE_H
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <tuple>
#include "protospacer_detection/array_context.h"
#include "protospacer_detection/spacer_boundaries.h"

namespace mcaat::proto {
using Node=uint64_t;
using Path=std::vector<Node>;
constexpr Node missing=std::numeric_limits<Node>::max();
inline std::string reverse_complement(std::string s) {
    std::reverse(s.begin(),s.end());
    for(char& c:s) {switch(c){case 'A':c='T';break;case 'T':c='A';break;case 'G':c='C';break;case 'C':c='G';break;default:throw std::runtime_error("non-ACGT sequence");}}
    return s;
}
inline std::string canonical(const std::string& s) {return std::min(s,reverse_complement(s));}
inline double median(std::vector<double> x) {
    if(x.empty())return 0;
    std::sort(x.begin(),x.end());return x.size()%2?x[x.size()/2]:(x[x.size()/2-1]+x[x.size()/2])/2;
}
struct AdapterPattern {
    std::string sequence;
    size_t exact_begin=0,exact_end=0,max_edits=0;
    bool allow_indels=false;
};
struct Options {
    size_t max_bp=1550,max_flank_bp=750,beam_width=32,max_states=50000,max_paths=16,max_seeds=128,min_match=0;
    size_t coverage_window=7,max_query_states=200000;
    size_t max_node_visits=3;
    size_t min_external_flank=0; // 0: k+1, one complete edge outside the seed
    bool joint_flank_check=true;
    bool resolve_graph_ties=false;
    bool junction_check=true,branch_guard=false;
    size_t branch_lookahead=75;
    double branch_margin=1.,prediction_margin=.1;
    BoundaryOptions boundaries;
    double max_coverage_penalty=20,read_weight=1,pair_weight=2,pam_weight=1,array_repeat_identity=.85;
    std::string coverage_model="auto";
    // Require the exact 13-base Illumina core AND confirming adapter sequence
    // (six model bases, <=2 edits). The short core alone can occur biologically.
    // Nextera and user-supplied adapters are exact. Empty explicitly disables.
    std::vector<AdapterPattern> adapters={{"AGATCGGAAGAGCACACGT",0,13,2,true},
                                         {"AGATCGGAAGAGCGTCGTG",0,13,2,true},
                                         {"CTGTCTCTTATACACATCT"}};
    void validate() const {
        boundaries.validate();
        if(!max_node_visits || max_node_visits>64)
            throw std::invalid_argument("max-node-visits must be in [1, 64]");
        if(!max_bp || !max_flank_bp || !branch_lookahead || branch_lookahead>1000 || !beam_width || !max_states || !max_paths || !max_seeds || !max_query_states || coverage_window<3)
            throw std::invalid_argument("extension sizes must be positive; coverage-window must be >=3");
        if(beam_width>1024 || max_paths>1024 || max_seeds>65536)
            throw std::invalid_argument("beam-width/max-paths must be <=1024; max-seeds <=65536");
        for(double d:{max_coverage_penalty,read_weight,pair_weight,pam_weight,branch_margin,prediction_margin})
            if(!std::isfinite(d) || d<0)throw std::invalid_argument("extension score settings must be finite and nonnegative");
        if(coverage_model!="auto" && coverage_model!="poisson" && coverage_model!="negative-binomial")
            throw std::invalid_argument("coverage-model must be auto, poisson or negative-binomial");
        if(!std::isfinite(array_repeat_identity) || array_repeat_identity<.5 || array_repeat_identity>1)
            throw std::invalid_argument("array-repeat-identity must be in [0.5, 1]");
        if(adapters.size()>64)throw std::invalid_argument("at most 64 adapter sequences are supported");
        for(const auto& p:adapters) {
            const auto& s=p.sequence;
            if(s.size()<12 || s.size()>256 || s.find_first_not_of("ACGT")!=std::string::npos)
                throw std::invalid_argument("adapter sequences must contain 12-256 ACGT bases");
            if(p.exact_begin>p.exact_end || p.exact_end>s.size() || p.max_edits>s.size()-(p.exact_end-p.exact_begin))
                throw std::invalid_argument("invalid adapter anchor or edit allowance");
            if(p.allow_indels && (p.max_edits>4 || p.exact_begin==p.exact_end || (p.exact_begin!=0 && p.exact_end!=s.size())))
                throw std::invalid_argument("adapter indels require a terminal exact anchor and at most four edits");
        }
    }
};
class AdapterGuard {
    std::vector<AdapterPattern> motifs_;
    size_t width_=0;
    static bool matches(const std::string& s,size_t start,const AdapterPattern& p) {
        size_t errors=0;
        for(size_t i=0;i<p.sequence.size();++i)if(s[start+i]!=p.sequence[i]) {
            if(i>=p.exact_begin && i<p.exact_end)return false;
            if(++errors>p.max_edits)return false;
        }
        return true;
    }
    static size_t edit_distance(const std::string& a,const std::string& b) {
        std::vector<size_t> row(b.size()+1),next(row.size());
        for(size_t j=0;j<row.size();++j)row[j]=j;
        for(size_t i=0;i<a.size();++i) {
            next[0]=i+1;
            for(size_t j=1;j<row.size();++j)next[j]=std::min({next[j-1]+1,row[j]+1,row[j-1]+(a[i]!=b[j-1])});
            row.swap(next);
        }
        return row.back();
    }
    static bool matches_span(const std::string& s,size_t start,size_t length,const AdapterPattern& p) {
        if(!p.allow_indels)return length==p.sequence.size() && matches(s,start,p);
        size_t anchor=p.exact_end-p.exact_begin;
        if(length<anchor)return false;
        if(p.exact_begin==0) {
            if(s.compare(start,anchor,p.sequence,0,anchor))return false;
            return edit_distance(p.sequence.substr(anchor),s.substr(start+anchor,length-anchor))<=p.max_edits;
        }
        if(s.compare(start+length-anchor,anchor,p.sequence,p.exact_begin,anchor))return false;
        return edit_distance(p.sequence.substr(0,p.exact_begin),s.substr(start,length-anchor))<=p.max_edits;
    }
public:
    AdapterGuard(const std::vector<AdapterPattern>& sequences,bool backward=false) {
        for(const auto& p:sequences)for(bool rc:{false,true}) {
            auto pattern=p;
            if(rc)pattern.sequence=reverse_complement(pattern.sequence);
            if(backward)std::reverse(pattern.sequence.begin(),pattern.sequence.end());
            if(rc!=backward) {
                pattern.exact_begin=p.sequence.size()-p.exact_end;
                pattern.exact_end=p.sequence.size()-p.exact_begin;
            }
            width_=std::max(width_,pattern.sequence.size()+(pattern.allow_indels?pattern.max_edits:0));motifs_.push_back(std::move(pattern));
        }
    }
    bool empty() const {return motifs_.empty();}
    bool contains(const std::string& s) const {
        for(const auto& m:motifs_) {
            size_t extra=m.allow_indels?m.max_edits:0;
            for(size_t n=m.sequence.size()-extra;n<=m.sequence.size()+extra;++n)
                for(size_t i=0;i+n<=s.size();++i)if(matches_span(s,i,n,m))return true;
        }
        return false;
    }
    std::string suffix(const std::string& s) const {return s.substr(s.size()>width_?s.size()-width_:0);}
    std::string append(std::string& tail,char base) const {
        if(empty())return {};
        tail+=base;if(tail.size()>width_)tail.erase(0,tail.size()-width_);
        std::string match;
        for(const auto& m:motifs_) {
            size_t extra=m.allow_indels?m.max_edits:0;
            for(size_t n=m.sequence.size()-extra;n<=m.sequence.size()+extra;++n)
                if(n<=tail.size() && n>match.size() && matches_span(tail,tail.size()-n,n,m))match=tail.substr(tail.size()-n);
        }
        return match;
    }
};
struct Arc {Node to,edge;double count;};
struct Seed {
    Path path;
    size_t begin=0,end=0; // half-open coordinates in the submitted spacer
    std::vector<Arc> left,right;
};
struct SeedSearch {std::vector<Seed> seeds;bool truncated=false;size_t eligible=0;};

// Keep distinct intervals and branches, even when they overlap other spacers.
// Optional edge exclusions are supported. Production source-array exclusion
// evaluates the complete local sequence context, so shared prefixes remain
// traversable until divergence or a sufficiently similar repeat is established.
template<class Graph>
SeedSearch isolate(Graph& graph,const Path& full,const std::set<Node>& array_edges,const Options& options) {
    SeedSearch result;if(full.empty())return result;
    std::unordered_set<Node> in_spacer(full.begin(),full.end());
    std::vector<std::vector<Arc>> left(full.size()),right(full.size());
    for(size_t i=0;i<full.size();++i) {
        for(auto arc:graph.neighbors(full[i],true))
            if(!in_spacer.count(arc.to) && !array_edges.count(arc.edge))left[i].push_back(arc);
        for(auto arc:graph.neighbors(full[i],false))
            if(!in_spacer.count(arc.to) && !array_edges.count(arc.edge))right[i].push_back(arc);
    }
    // Longer exact matches first. One-sided intervals must reach the corresponding
    // original spacer end; internal fragments require external evidence on both ends.
    for(size_t span=full.size();span;--span)for(size_t i=0;i+span<=full.size();++i) {
        size_t j=i+span-1,matched=span+graph.k()-1;
        if(matched<std::max(graph.k(),options.min_match))continue;
        if((left[i].empty() && i!=0) || (right[j].empty() && j+1!=full.size()))continue;
        if(left[i].empty() && right[j].empty())continue;
        ++result.eligible;
        if(result.seeds.size()>=options.max_seeds){result.truncated=true;continue;}
        result.seeds.push_back({Path(full.begin()+i,full.begin()+j+1),i,j+graph.k(),left[i],right[j]});
    }
    return result;
}

struct MultiplicityProfile {
    size_t samples=0,weak=0;
    double mean_log=0,m2=0,step_sum=0,previous=0;
    void observe(double count) {
        if(count<=1)++weak;
        double value=std::log2(1+std::max(0.,count));
        if(samples)step_sum+=(value-previous)*(value-previous);
        ++samples;double delta=value-mean_log;mean_log+=delta/samples;m2+=delta*(value-mean_log);previous=value;
    }
    double variance() const {return samples?std::max(0.,m2/samples):0;}
    double roughness() const {return samples>1?step_sum/(samples-1):0;}
    double level() const {return samples?std::exp2(mean_log)-1:0;}
    double weak_fraction() const {return samples?double(weak)/samples:0;}
    double loss() const {return variance()+roughness()+weak_fraction()+1./(samples+1);}
};
struct Coverage {
    std::vector<double> baseline;
    std::vector<double> recent;
    MultiplicityProfile profile;
    double penalty=0;size_t scored=0;
    double mean() const {return median(baseline);}
    static double dispersion_of(const std::vector<double>& values) {
        if(values.size()<3)return 0;
        const double m=median(values);std::vector<double> delta;
        for(double c:values)delta.push_back(std::abs(c-m));
        double sigma=1.4826*median(delta);
        return m>0?std::max(0.,(sigma*sigma-m)/(m*m)):0;
    }
    double dispersion() const {return dispersion_of(baseline);}
    double local_mean() const {return recent.empty()?mean():median(recent);}
    double local_dispersion() const {return recent.empty()?dispersion():dispersion_of(recent);}
    // Log-likelihood loss relative to a saturated mean, evaluated at separated
    // count observations. This is a ranking score, not a calibrated posterior.
    static double loss(double count,double mean,double dispersion) {
        mean=std::max(mean,1e-6);count=std::max(count,0.);
        if(dispersion<=1e-9)return std::max(0.,mean-count+(count?count*std::log(count/mean):0));
        double r=1/dispersion;
        return std::max(0.,(count?count*std::log(count/mean):0)-(count+r)*std::log((count+r)/(mean+r)));
    }
    double observe(double count,const Options& options,bool sample=true) {
        profile.observe(count);
        if(!sample)return 0;
        if(baseline.size()<options.coverage_window) {baseline.push_back(count);recent.push_back(count);return 0;}
        double alpha=options.coverage_model=="poisson"?0:local_dispersion();
        // Explicit NB uses a weak overdispersion floor when local MAD is zero.
        if(options.coverage_model=="negative-binomial")alpha=std::max(alpha,.05);
        double value=loss(count,local_mean(),alpha);penalty+=value;++scored;
        // Test before adapting: a rejected jump cannot move the local baseline.
        // Keep the initial baseline separately for the two-flank comparison.
        if(value<=options.max_coverage_penalty) {
            if(recent.size()>=options.coverage_window)recent.erase(recent.begin());
            recent.push_back(count);
        }
        return value;
    }
    double score() const {return scored?-penalty/scored:0;}
};
struct Stability {
    double within=0,jumps=0,between=0,support=0,weak=0;
    double score() const {return -(within+jumps+between+support+weak);}
};
inline Stability stability(const Coverage& left,const Coverage& right) {
    const auto& l=left.profile;const auto& r=right.profile;
    Stability s;
    s.within=(l.variance()+r.variance())/2;
    s.jumps=(l.roughness()+r.roughness())/2;
    s.weak=(l.weak_fraction()+r.weak_fraction())/2;
    if(l.samples && r.samples)s.between=(l.mean_log-r.mean_log)*(l.mean_log-r.mean_log);
    // A short/absent flank must not earn a perfect score from zero observations.
    s.support=1./(l.samples+1)+1./(r.samples+1);return s;
}
struct Side {
    Path nodes; // starts at seed endpoint, grows outwards
    Coverage coverage;
    MultiplicityProfile beam_profile; // provisional: includes unresolved boundary edges
    size_t source_run=0;
    size_t branches=0,since_sample=0;
    std::string stop;
    std::string context="unexamined",context_prefix;
    size_t context_clear_at=0;
    std::string adapter_tail;
    bool adapter_overlaps_seed=false;
};
struct AdapterStop {
    bool backward=false,overlaps_seed=false;
    std::string motif; // in the query's forward orientation, even on the left
    size_t clean_flank_bp=0;
    Node clean_endpoint=missing,detected_node=missing;
};
struct BranchCheck {
    bool backward=false,chosen=false;
    Node vertex=missing,edge=missing,target=missing;
    size_t observations=0,steps=0;
    double score=0;
    std::string status,decision;
};
struct SearchDiagnostics {
    size_t adapter_stops=0;
    std::vector<AdapterStop> adapters;
    size_t branch_checks=0;
    std::vector<BranchCheck> branches;
    // Bound diagnostics independently of search work. Counts remain complete.
    void record(AdapterStop stop) {++adapter_stops;if(adapters.size()<4096)adapters.push_back(std::move(stop));}
    void record_branch(BranchCheck check) {++branch_checks;if(branches.size()<4096)branches.push_back(std::move(check));}
};
template<class Graph>
Coverage path_coverage(Graph& graph,const Path& nodes,size_t length,bool backwards,size_t clear_at,
                       const std::set<Node>* shared_context,const Options& options) {
    Coverage cov;size_t since=0;
    for(size_t i=1;i<=length;++i) {
        if(i<clear_at)continue;
        Arc chosen{missing,missing,0};
        for(auto a:graph.neighbors(nodes[i-1],backwards))if(a.to==nodes[i]){chosen=a;break;}
        if(chosen.to==missing)throw std::runtime_error("missing coverage path edge");
        if(shared_context && shared_context->count(chosen.edge)){++since;continue;}
        bool junction=options.junction_check &&
            (graph.neighbors(nodes[i-1],backwards).size()>1 || graph.neighbors(nodes[i],!backwards).size()>1);
        bool sample=junction || cov.baseline.empty() || since>=graph.k()+1;
        cov.observe(chosen.count,options,sample);
        if(sample)since=0;else ++since;
    }
    return cov;
}
struct WalkResult {std::vector<Side> paths;size_t expanded=0;bool beam_pruned=false,work_limited=false;};
inline bool better_side(const Side& a,const Side& b) {
    if(a.adapter_overlaps_seed!=b.adapter_overlaps_seed)return !a.adapter_overlaps_seed;
    if((a.nodes.size()>1)!=(b.nodes.size()>1))return a.nodes.size()>1;
    if(a.coverage.profile.loss()!=b.coverage.profile.loss())return a.coverage.profile.loss()<b.coverage.profile.loss();
    if(a.nodes.size()!=b.nodes.size())return a.nodes.size()>b.nodes.size();
    return a.nodes<b.nodes;
}

struct BranchProbe {
    bool viable=true,limited=false;
    size_t observations=0,steps=0;
    double score=0;
    std::string status="insufficient_evidence";
};
// A shared vertex alone is not evidence of an array return. Reject a closed
// walk only when every edge since the previous visit belongs to source cycles.
// Other repeats are bounded by visits per vertex and the ordinary work/bp caps.
class SourceCycleGuard {
    std::unordered_set<std::string> words_;
    std::unordered_map<Node,bool> cache_;
public:
    SourceCycleGuard(size_t k,const std::vector<FlankContext>* models) {
        if(models)for(const auto& m:*models)if(!m.cycle.empty()) {
            std::string wrapped=m.cycle;
            while(wrapped.size()<m.cycle.size()+k)wrapped+=m.cycle;
            for(size_t i=0;i<m.cycle.size();++i)words_.insert(wrapped.substr(i,k+1));
        }
    }
    template<class Graph> bool contains(Graph& graph,Node from,const Arc& arc,bool back) {
        if(words_.empty())return false;
        auto it=cache_.find(arc.edge);if(it!=cache_.end())return it->second;
        auto label=graph.label(from);char base=graph.extension_base(arc,back);
        bool found=words_.count(back?std::string(1,base)+label:label+base);
        cache_.emplace(arc.edge,found);return found;
    }
    template<class Graph> bool contains(Graph& graph,Node from,Node to,bool back) {
        for(auto arc:graph.neighbors(from,back))if(arc.to==to)return contains(graph,from,arc,back);
        throw std::runtime_error("missing source-cycle path edge");
    }
};
inline std::string repeat_stop(const Path& seed,const Side& state,Node next,size_t source_run,size_t limit) {
    size_t count=0,seed_count=0,last=0;
    for(size_t i=0;i<seed.size();++i)if(seed[i]==next){++count;++seed_count;last=i;}
    for(size_t i=1;i<state.nodes.size();++i)if(state.nodes[i]==next){++count;last=seed.size()-1+i;}
    if(!count)return {};
    size_t position=seed.size()-1+state.nodes.size();
    if(source_run>=position-last)return "source_array_cycle";
    if(count>=std::max(size_t(1),seed_count)+limit-1)return "repeat_visit_limit";
    return {};
}
template<class Graph>
bool compatible_walk(Graph& graph,const Path& path,const Path& seed,size_t seed_begin,
                     SourceCycleGuard& guard,const Options& options) {
    std::unordered_map<Node,size_t> originals,counts,last;
    for(auto v:seed)++originals[v];
    size_t run=0,seed_end=seed_begin+seed.size();
    for(size_t i=0;i<path.size();++i) {
        if(i)run=guard.contains(graph,path[i-1],path[i],false)?run+1:0;
        auto v=path[i];
        if(++counts[v]>std::max(size_t(1),originals[v])+options.max_node_visits-1)return false;
        auto previous=last.find(v);
        if(previous!=last.end() && run>=i-previous->second &&
           !(previous->second>=seed_begin && i<seed_end))return false;
        last[v]=i;
    }
    return true;
}
template<class Graph>
BranchProbe probe_branch(Graph& graph,const Side& current,Arc first,bool backward,
                         const Path& oriented_seed,SourceCycleGuard& source_guard,const std::set<Node>& masked,
                         const std::set<Node>* shared_context,const std::vector<FlankContext>* contexts,
                         const AdapterGuard& guard,size_t budget,size_t& expanded,
                         const Options& options,SearchDiagnostics* diagnostics) {
    BranchProbe result;Side walk=current;Path& path=walk.nodes;
    std::string prefix=current.context_prefix,context=current.context,tail=current.adapter_tail;
    MultiplicityProfile local;double loss=0;
    const bool baseline=current.coverage.baseline.size()>=3;
    double alpha=options.coverage_model=="poisson"?0:current.coverage.local_dispersion();
    if(options.coverage_model=="negative-binomial")alpha=std::max(alpha,.05);
    Arc arc=first;
    for(size_t step=0;step<std::min(budget,options.branch_lookahead);++step) {
        if(expanded>=options.max_states){result.limited=true;result.status="work_limit";return result;}
        ++expanded;++result.steps;
        size_t source_run=source_guard.contains(graph,path.back(),arc,backward)?walk.source_run+1:0;
        auto stop=repeat_stop(oriented_seed,walk,arc.to,source_run,options.max_node_visits);
        if(!stop.empty() || masked.count(arc.edge)){result.status=stop.empty()?"source_array_context":stop;break;}
        char base=graph.extension_base(arc,backward);
        auto adapter=guard.append(tail,base);
        if(!adapter.empty()) {
            size_t extended=path.size(),keep=extended>adapter.size()?extended-adapter.size():0;
            if(backward)std::reverse(adapter.begin(),adapter.end());
            if(diagnostics)diagnostics->record({backward,adapter.size()>extended,adapter,keep,path[keep],arc.to});
            result.viable=false;result.status="sequencing_adapter";return result;
        }
        if(context=="pending") {
            prefix+=base;auto state=context_state(prefix,*contexts,options.array_repeat_identity);
            if(state==ContextState::array){result.viable=false;result.status="source_array_context";return result;}
            if(state==ContextState::external)context="external";
        }
        if(context!="pending" && (!shared_context || !shared_context->count(arc.edge))) {
            double penalty=baseline?Coverage::loss(arc.count,current.coverage.local_mean(),alpha):0;
            if(baseline && penalty>options.max_coverage_penalty) {
                result.viable=false;result.status="coverage_break";return result;
            }
            loss+=penalty;local.observe(arc.count);++result.observations;
        }
        path.push_back(arc.to);walk.source_run=source_run;
        auto next=graph.neighbors(arc.to,backward);
        bool reaches_mask=std::any_of(next.begin(),next.end(),[&](const auto& a){return masked.count(a.edge);});
        next.erase(std::remove_if(next.begin(),next.end(),[&](const auto& a){
            size_t run=source_guard.contains(graph,path.back(),a,backward)?walk.source_run+1:0;
            return !repeat_stop(oriented_seed,walk,a.to,run,options.max_node_visits).empty()||masked.count(a.edge);
        }),next.end());
        if(next.empty() && reaches_mask){result.viable=false;result.status="source_array_context";return result;}
        if(next.size()!=1)break;
        arc=next.front();
    }
    if(context=="pending")result.status="incomplete_context";
    else if(result.observations>=3) {
        // A mean loss, not summed support: more correlated k-mers must not
        // manufacture confidence or reward a longer branch.
        result.score=loss/result.observations+local.variance()+local.roughness();
        result.status="scored";
    }
    return result;
}

template<class Graph>
WalkResult extend(Graph& graph,const Seed& seed,bool backward,size_t budget,
                  const std::set<Node>& array_edges,const Options& options,const std::set<Node>* shared_context=nullptr,
                  const std::vector<FlankContext>* contexts=nullptr,SearchDiagnostics* diagnostics=nullptr) {
    WalkResult result;
    const auto& first=backward?seed.left:seed.right;
    Node endpoint=backward?seed.path.front():seed.path.back();
    Side initial;initial.nodes={endpoint};
    Path oriented_seed=seed.path;if(backward)std::reverse(oriented_seed.begin(),oriented_seed.end());
    SourceCycleGuard source_guard(graph.k(),contexts);
    for(size_t i=1;i<oriented_seed.size();++i)
        initial.source_run=source_guard.contains(graph,oriented_seed[i-1],oriented_seed[i],backward)?initial.source_run+1:0;
    AdapterGuard adapter_guard(options.adapters,backward);
    if(!adapter_guard.empty()) {
        auto sequence=graph.spell(seed.path);if(backward)std::reverse(sequence.begin(),sequence.end());
        initial.adapter_tail=adapter_guard.suffix(sequence);
    }
    if(contexts && !contexts->empty())initial.context="pending";
    if(!budget || first.empty()) {
        initial.stop=budget?"no_external_flank":"length_limit";
        if(!budget && !first.empty())initial.context="incomplete";
        result.paths={initial};return result;
    }
    std::vector<Side> frontier{initial};
    auto finish=[&](Side state,std::string reason) {
        if(state.context=="pending")state.context=reason=="array_context"?"array":"incomplete";
        if(state.context=="array" || state.context=="incomplete") {
            state.nodes.resize(1);state.coverage={};
            if(state.context=="array")reason="source_array_context";
        }
        state.stop=std::move(reason);
        for(const auto& old:result.paths)if(old.nodes==state.nodes && old.stop==state.stop && old.context==state.context && old.adapter_overlaps_seed==state.adapter_overlaps_seed)return;
        result.paths.push_back(std::move(state));
        // Bound terminal storage as well as the frontier.
        if(result.paths.size()>options.beam_width*2) {
            std::sort(result.paths.begin(),result.paths.end(),better_side);
            result.paths.resize(options.beam_width);result.beam_pruned=true;
        }
    };
    // Level-synchronous beam: every live path receives one more base before
    // pruning. Terminal paths never displace still-growing alternatives.
    while(!frontier.empty()) {
        std::vector<Side> next_frontier;
        for(auto& entry:frontier) {
            Side current=std::move(entry);
            if(current.nodes.size()-1>=budget){finish(std::move(current),"length_limit");continue;}
            if(result.expanded>=options.max_states) {
                result.work_limited=true;finish(std::move(current),"work_limit");continue;
            }
            auto arcs=current.nodes.size()==1?first:graph.neighbors(current.nodes.back(),backward);
            std::vector<Arc> usable;bool masked=false;std::string repeat_reason;
            for(auto a:arcs) {
                size_t run=source_guard.contains(graph,current.nodes.back(),a,backward)?current.source_run+1:0;
                auto stop=repeat_stop(oriented_seed,current,a.to,run,options.max_node_visits);
                if(!stop.empty()){repeat_reason=stop;continue;}
                if(array_edges.count(a.edge)){masked=true;continue;}
                usable.push_back(a);
            }
            std::sort(usable.begin(),usable.end(),[](auto a,auto b){return a.to<b.to;});
            if(usable.empty()){finish(std::move(current),!repeat_reason.empty()?repeat_reason:masked?"array_context":"dead_end");continue;}
            bool branch=usable.size()>1;
            if(branch)++current.branches;
            if(branch && options.branch_guard) {
                std::vector<BranchProbe> probes;std::vector<size_t> viable;
                bool limited=false;
                for(size_t j=0;j<usable.size();++j) {
                    probes.push_back(probe_branch(graph,current,usable[j],backward,oriented_seed,source_guard,array_edges,shared_context,
                        contexts,adapter_guard,budget-(current.nodes.size()-1),result.expanded,options,diagnostics));
                    limited|=probes.back().limited;
                    if(probes.back().viable)viable.push_back(j);
                }
                size_t chosen=usable.size();std::string decision="ambiguous_branch";
                if(limited){result.work_limited=true;decision="work_limit";}
                else if(viable.empty())decision=std::all_of(probes.begin(),probes.end(),[](const auto& p){
                    return p.status=="source_array_context";
                })?"array_context":"no_supported_branch";
                else if(viable.size()==1) {chosen=viable.front();decision="only_supported_branch";}
                else {
                    bool scored=std::all_of(viable.begin(),viable.end(),[&](auto j){return probes[j].status=="scored";});
                    std::sort(viable.begin(),viable.end(),[&](auto a,auto b){
                        return probes[a].score!=probes[b].score?probes[a].score<probes[b].score:usable[a].to<usable[b].to;
                    });
                    if(scored && probes[viable[1]].score-probes[viable[0]].score>options.branch_margin) {
                        chosen=viable[0];decision="lookahead_winner";
                    }
                }
                if(diagnostics)for(size_t j=0;j<usable.size();++j)
                    diagnostics->record_branch({backward,j==chosen,current.nodes.back(),usable[j].edge,usable[j].to,
                        probes[j].observations,probes[j].steps,probes[j].score,probes[j].status,decision});
                if(chosen==usable.size()) {finish(std::move(current),decision);continue;}
                usable={usable[chosen]};
            }
            std::vector<Side> children;
            std::vector<std::pair<Side,std::string>> rejected;
            for(const auto& a:usable) {
                if(result.expanded>=options.max_states){result.work_limited=true;break;}
                ++result.expanded;
                Side next=branch?current:std::move(current);
                next.source_run=source_guard.contains(graph,next.nodes.back(),a,backward)?next.source_run+1:0;
                next.beam_profile.observe(a.count);
                char base=0;
                if(!adapter_guard.empty() || next.context=="pending")base=graph.extension_base(a,backward);
                auto adapter=adapter_guard.append(next.adapter_tail,base);
                if(!adapter.empty()) {
                    // The completing base has not been appended to nodes yet.
                    // Roll back the entire motif, not just its final nucleotide.
                    const size_t extended=next.nodes.size();
                    next.adapter_overlaps_seed=adapter.size()>extended;
                    size_t keep=extended>adapter.size()?extended-adapter.size():0;
                    next.nodes.resize(keep+1);
                    if(backward)std::reverse(adapter.begin(),adapter.end());
                    if(diagnostics)diagnostics->record({backward,next.adapter_overlaps_seed,adapter,keep,next.nodes.back(),a.to});
                    if(next.context_prefix.size()>keep)next.context_prefix.resize(keep);
                    if(contexts && !contexts->empty() && keep<next.context_clear_at) {
                        auto state=context_state(next.context_prefix,*contexts,options.array_repeat_identity);
                        next.context=state==ContextState::external?"external":state==ContextState::array?"array":"pending";
                    }
                    next.coverage=path_coverage(graph,next.nodes,keep,backward,next.context_clear_at,shared_context,options);
                    rejected.emplace_back(std::move(next),"sequencing_adapter");continue;
                }
                if(next.context=="pending") {
                    next.context_prefix+=base;
                    auto state=context_state(next.context_prefix,*contexts,options.array_repeat_identity);
                    if(state==ContextState::array) {
                        next.context="array";rejected.emplace_back(std::move(next),"source_array_context");continue;
                    }
                    if(state==ContextState::external) {
                        next.context="external";next.context_clear_at=next.nodes.size();
                    }
                }
                bool junction=options.junction_check && (arcs.size()>1 || graph.neighbors(a.to,!backward).size()>1);
                bool sample=(!shared_context || !shared_context->count(a.edge)) &&
                    next.context!="pending" &&
                    (junction || next.coverage.baseline.empty() || next.since_sample>=graph.k()+1);
                auto before=next.coverage;
                bool evidence=next.context!="pending" && (!shared_context || !shared_context->count(a.edge));
                double penalty=evidence?next.coverage.observe(a.count,options,sample):0;
                if(sample)next.since_sample=0;else ++next.since_sample;
                if(penalty>options.max_coverage_penalty) {next.coverage=std::move(before);rejected.emplace_back(std::move(next),"coverage_break");continue;}
                next.nodes.push_back(a.to);children.push_back(std::move(next));
            }
            // Rejected siblings cannot manufacture a short winning prefix.
            if(children.empty()) {
                if(rejected.empty())finish(std::move(current),"work_limit");
                else for(auto& stop:rejected)finish(std::move(stop.first),std::move(stop.second));
            }
            for(auto& child:children)next_frontier.push_back(std::move(child));
        }
        if(next_frontier.size()>options.beam_width) {
            std::sort(next_frontier.begin(),next_frontier.end(),[](const Side& a,const Side& b){
                if(a.beam_profile.loss()!=b.beam_profile.loss())return a.beam_profile.loss()<b.beam_profile.loss();
                return a.nodes<b.nodes;
            });
            next_frontier.resize(options.beam_width);result.beam_pruned=true;
        }
        frontier=std::move(next_frontier);
    }
    std::sort(result.paths.begin(),result.paths.end(),better_side);
    // The same partial path can stop before several incompatible successors.
    std::set<std::tuple<Path,std::string,bool>> seen;std::vector<Side> unique;
    for(auto& side:result.paths)if(seen.emplace(side.nodes,side.context,side.adapter_overlaps_seed).second)unique.push_back(std::move(side));
    if(unique.size()>options.beam_width){unique.resize(options.beam_width);result.beam_pruned=true;}
    result.paths=std::move(unique);return result;
}
struct Candidate {
    Path path;
    Side left,right;
    size_t seed_begin=0,seed_end=0,expanded=0,left_length=0,right_length=0,external_left=0,external_right=0;
    bool search_limited=false;
    double score=0,pam_bonus=0;
    std::string sequence;
    size_t junctions=0,supported_junctions=0,read_support=0,pair_support=0,unique_pairs=0;
    bool short_flank_trimmed=false,flank_conflict=false;
    bool ambiguity_trimmed=false;
    double joint_flank_loss=0;
};
struct CommonFlanks {size_t left=0,right=0;bool seed_agrees=true;};
inline CommonFlanks common_flanks(const Candidate& winner,size_t winner_query_begin,
                                 const Candidate& alternative,size_t alternative_query_begin) {
    CommonFlanks result;
    const int64_t aligned=int64_t(alternative.seed_begin)+int64_t(winner_query_begin)-int64_t(alternative_query_begin);
    size_t seed_length=winner.seed_end-winner.seed_begin;
    if(aligned<0 || uint64_t(aligned)+seed_length>alternative.sequence.size() ||
       alternative.sequence.compare(size_t(aligned),seed_length,winner.sequence,winner.seed_begin,seed_length)) {
        result.seed_agrees=false;return result;
    }
    size_t begin=size_t(aligned),end=begin+seed_length;
    while(result.left<winner.seed_begin && result.left<begin &&
          winner.sequence[winner.seed_begin-result.left-1]==alternative.sequence[begin-result.left-1])++result.left;
    while(winner.seed_end+result.right<winner.sequence.size() && end+result.right<alternative.sequence.size() &&
          winner.sequence[winner.seed_end+result.right]==alternative.sequence[end+result.right])++result.right;
    return result;
}
inline void retain_flanks(Candidate& c,size_t left,size_t right,size_t k) {
    if(left>c.seed_begin || right>c.sequence.size()-c.seed_end)throw std::runtime_error("invalid retained flank length");
    size_t begin=c.seed_begin-left,end=c.seed_end+right;
    c.path=Path(c.path.begin()+begin,c.path.begin()+end-k+1);
    c.sequence=c.sequence.substr(begin,end-begin);c.seed_begin-=begin;c.seed_end-=begin;
    c.left_length=left;c.right_length=right;
}
inline void trim_left(Candidate& c) {
    size_t remove=c.seed_begin;c.path.erase(c.path.begin(),c.path.begin()+remove);
    c.sequence.erase(0,remove);c.seed_end-=remove;c.seed_begin=0;c.left_length=0;
    c.left.coverage={};c.external_left=0;
}
inline void trim_right(Candidate& c,size_t k) {
    c.sequence.resize(c.seed_end);c.path.resize(c.seed_end-k+1);c.right_length=0;
    c.right.coverage={};c.external_right=0;
}
inline double joint_flank_loss(const Coverage& left,const Coverage& right,const Options& options) {
    // Insufficient observations cannot establish a conflict. Do not compare
    // summed spacer depth against flank depth: source and target share the seed.
    if(left.baseline.size()<3 || right.baseline.size()<3)return 0;
    double alpha=options.coverage_model=="poisson"?0:std::max(left.dispersion(),right.dispersion());
    if(options.coverage_model=="negative-binomial")alpha=std::max(alpha,.05);
    return std::max(Coverage::loss(left.mean(),right.mean(),alpha),Coverage::loss(right.mean(),left.mean(),alpha));
}
inline std::vector<Candidate> refine_flanks(Candidate c,size_t k,const Options& options) {
    size_t minimum=options.min_external_flank?options.min_external_flank:k+1;
    if(c.left_length && c.left_length<minimum) {
        trim_left(c);c.left.stop="insufficient_external_sequence";c.left.context="incomplete";c.short_flank_trimmed=true;
    }
    if(c.right_length && c.right_length<minimum) {
        trim_right(c,k);c.right.stop="insufficient_external_sequence";c.right.context="incomplete";c.short_flank_trimmed=true;
    }
    c.joint_flank_loss=joint_flank_loss(c.left.coverage,c.right.coverage,options);
    if(options.joint_flank_check && c.joint_flank_loss>options.max_coverage_penalty) {
        // Both sides may be real contexts at different abundances. Retain them
        // as separate candidates; never force them into one reconstructed locus.
        c.flank_conflict=true;Candidate right=c;
        trim_right(c,k);c.right.stop="flank_multiplicity_conflict";c.right.context="not_retained";
        trim_left(right);right.left.stop="flank_multiplicity_conflict";right.left.context="not_retained";
        c.score=stability(c.left.coverage,c.right.coverage).score();
        right.score=stability(right.left.coverage,right.right.coverage).score();
        return {std::move(c),std::move(right)};
    }
    c.score=stability(c.left.coverage,c.right.coverage).score();return {std::move(c)};
}
inline bool better_prediction(const Candidate& a,size_t a_matched,const Candidate& b,size_t b_matched,const Options& options) {
    if(options.pair_weight>0 && bool(a.pair_support)!=bool(b.pair_support))return bool(a.pair_support);
    if(a.score!=b.score)return a.score>b.score;
    if(a_matched!=b_matched)return a_matched>b_matched;
    auto an=a.left.coverage.profile.samples+a.right.coverage.profile.samples;
    auto bn=b.left.coverage.profile.samples+b.right.coverage.profile.samples;
    if(an!=bn)return an>bn;
    return a.sequence.size()!=b.sequence.size()?a.sequence.size()>b.sequence.size():a.sequence<b.sequence;
}
class RetainedIntervalContexts {
    struct Kept {size_t begin,end;std::string id;};
    std::map<std::pair<std::string,int64_t>,std::vector<Kept>> contexts_;
public:
    // Call in decreasing matched-interval length, AFTER output limits.
    std::string replacement(const std::string& sequence,size_t seed_begin,
                            size_t begin,size_t end,const std::string& id) {
        auto& previous=contexts_[{sequence,int64_t(seed_begin)-int64_t(begin)}];
        for(const auto& parent:previous)
            if(parent.begin<=begin && parent.end>=end)return parent.id;
        previous.push_back({begin,end,id});return {};
    }
};
template<class Graph>
std::vector<Candidate> reconstruct(Graph& graph,const Seed& seed,const std::set<Node>& array_edges,const Options& options,const std::set<Node>* shared_context=nullptr,
                                   const std::vector<FlankContext>* left_context=nullptr,const std::vector<FlankContext>* right_context=nullptr,
                                   SearchDiagnostics* diagnostics=nullptr) {
    size_t seed_bp=seed.path.size()+graph.k()-1;
    size_t budget=options.max_bp>seed_bp?options.max_bp-seed_bp:0;
    // Each side can use the whole budget; combination trims symmetrically and
    // transfers spare capacity, so a blocked side never erases the other flank.
    if(AdapterGuard(options.adapters).contains(graph.spell(seed.path))) {
        Candidate c;c.path=seed.path;c.sequence=graph.spell(seed.path);c.seed_end=seed_bp;
        c.left.stop=c.right.stop="adapter_in_seed";return {c};
    }
    auto side_budget=std::min(budget,options.max_flank_bp);
    auto left=extend(graph,seed,true,side_budget,array_edges,options,shared_context,left_context,diagnostics);
    auto right=extend(graph,seed,false,side_budget,array_edges,options,shared_context,right_context,diagnostics);
    std::vector<FlankContext> cycle_models;
    if(left_context)cycle_models.insert(cycle_models.end(),left_context->begin(),left_context->end());
    if(right_context)cycle_models.insert(cycle_models.end(),right_context->begin(),right_context->end());
    SourceCycleGuard join_guard(graph.k(),&cycle_models);
    std::vector<Candidate> results;std::set<std::tuple<Path,std::string,std::string>> seen;size_t combinations=0;
    bool combination_limit=left.paths.size()*right.paths.size()>options.beam_width*4;
    for(size_t diagonal=0;diagonal<left.paths.size()+right.paths.size();++diagonal)
      for(size_t li=0;li<left.paths.size() && li<=diagonal;++li) {
        size_t ri=diagonal-li;if(ri>=right.paths.size())continue;
        if(combinations++>=options.beam_width*4)continue;
        const auto& l=left.paths[li];const auto& r=right.paths[ri];
        if(l.adapter_overlaps_seed || r.adapter_overlaps_seed)continue;
        size_t ln=std::min(l.nodes.size()-1,budget/2);
        size_t rn=std::min(r.nodes.size()-1,budget-ln);
        ln=std::min(l.nodes.size()-1,budget-rn);
        Candidate c;c.left=l;c.right=r;
        if(ln+1<l.nodes.size())c.left.stop="length_limit";
        if(rn+1<r.nodes.size())c.right.stop="length_limit";
        c.left_length=ln;c.right_length=rn;
        // Recompute evidence only on emitted flanks, excluding trimmed suffixes.
        c.left.coverage=path_coverage(graph,l.nodes,ln,true,l.context_clear_at,shared_context,options);
        c.right.coverage=path_coverage(graph,r.nodes,rn,false,r.context_clear_at,shared_context,options);
        Path{}.swap(c.left.nodes);Path{}.swap(c.right.nodes);
        for(size_t i=ln;i>0;--i)c.path.push_back(l.nodes[i]);
        c.seed_begin=ln;c.seed_end=ln+seed_bp;
        c.path.insert(c.path.end(),seed.path.begin(),seed.path.end());
        for(size_t i=1;i<=rn;++i)c.path.push_back(r.nodes[i]);
        // The same bounds apply across the joined flanks. Existing seed repeats
        // remain intact; new source-array loops are never permitted.
        if(!compatible_walk(graph,c.path,seed.path,ln,join_guard,options))continue;
        if(!seen.emplace(c.path,c.left.context,c.right.context).second)continue;
        c.sequence=graph.spell(c.path);
        if(c.sequence.size()!=c.path.size()+graph.k()-1)throw std::runtime_error("incorrect reconstructed length");
        c.score=stability(c.left.coverage,c.right.coverage).score();
        c.expanded=left.expanded+right.expanded;
        c.search_limited=left.beam_pruned||right.beam_pruned||left.work_limited||right.work_limited||combination_limit;
        results.push_back(std::move(c));
    }
    if(results.size()>options.beam_width) {
        std::sort(results.begin(),results.end(),[&](const auto& a,const auto& b){
            if((a.sequence.size()>seed_bp)!=(b.sequence.size()>seed_bp))return a.sequence.size()>seed_bp;
            if(a.score!=b.score)return a.score>b.score;
            return a.sequence.size()!=b.sequence.size()?a.sequence.size()>b.sequence.size():a.sequence<b.sequence;
        });
        results.resize(options.beam_width);for(auto& c:results)c.search_limited=true;
    }
    // Never erase the seed if two otherwise valid flanks intersect.
    if(results.empty()) {
        Candidate c;c.path=seed.path;c.sequence=graph.spell(c.path);c.seed_end=seed_bp;
        c.left.stop=c.right.stop="incompatible_flanks";c.search_limited=true;
        auto contaminated=[](const WalkResult& walk){return !walk.paths.empty() && std::all_of(walk.paths.begin(),walk.paths.end(),[](const auto& s){return s.adapter_overlaps_seed;});};
        if(contaminated(left) || contaminated(right)) {
            c.left.stop=c.right.stop="adapter_in_seed";
            c.search_limited=left.beam_pruned||right.beam_pruned||left.work_limited||right.work_limited||combination_limit;
        }
        c.expanded=left.expanded+right.expanded;
        results.push_back(std::move(c));
    }
    return results;
}
} // namespace mcaat::proto
#endif
