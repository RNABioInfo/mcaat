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
#include <vector>
#include <tuple>
#include "protospacer_detection/array_context.h"

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
struct Options {
    size_t max_bp=2500,beam_width=32,max_states=50000,max_paths=16,max_seeds=128,min_match=0;
    size_t coverage_window=7,max_query_states=200000;
    double max_coverage_penalty=20,read_weight=1,pair_weight=2,pam_weight=1,array_repeat_identity=.85;
    std::string coverage_model="auto";
    void validate() const {
        if(!max_bp || !beam_width || !max_states || !max_paths || !max_seeds || !max_query_states || coverage_window<3)
            throw std::invalid_argument("extension sizes must be positive; coverage-window must be >=3");
        if(beam_width>1024 || max_paths>1024 || max_seeds>65536)
            throw std::invalid_argument("beam-width/max-paths must be <=1024; max-seeds <=65536");
        for(double d:{max_coverage_penalty,read_weight,pair_weight,pam_weight})
            if(!std::isfinite(d) || d<0)throw std::invalid_argument("extension score settings must be finite and nonnegative");
        if(coverage_model!="auto" && coverage_model!="poisson" && coverage_model!="negative-binomial")
            throw std::invalid_argument("coverage-model must be auto, poisson or negative-binomial");
        if(!std::isfinite(array_repeat_identity) || array_repeat_identity<.5 || array_repeat_identity>1)
            throw std::invalid_argument("array-repeat-identity must be in [0.5, 1]");
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
    size_t samples=0;
    double mean_log=0,m2=0,step_sum=0,previous=0;
    void observe(double count) {
        double value=std::log2(1+std::max(0.,count));
        if(samples)step_sum+=(value-previous)*(value-previous);
        ++samples;double delta=value-mean_log;mean_log+=delta/samples;m2+=delta*(value-mean_log);previous=value;
    }
    double variance() const {return samples?std::max(0.,m2/samples):0;}
    double roughness() const {return samples>1?step_sum/(samples-1):0;}
    double level() const {return samples?std::exp2(mean_log)-1:0;}
    double loss() const {return variance()+roughness()+1./(samples+1);}
};
struct Coverage {
    std::vector<double> baseline;
    MultiplicityProfile profile;
    double penalty=0;size_t scored=0;
    double mean() const {return median(baseline);}
    double dispersion() const {
        if(baseline.size()<3)return 0;
        const double m=mean();std::vector<double> delta;
        for(double c:baseline)delta.push_back(std::abs(c-m));
        double sigma=1.4826*median(delta);
        return m>0?std::max(0.,(sigma*sigma-m)/(m*m)):0;
    }
    // Log-likelihood loss relative to a saturated mean, evaluated at separated
    // count observations. This is a ranking score, not a calibrated posterior.
    static double loss(double count,double mean,double dispersion) {
        mean=std::max(mean,1e-6);count=std::max(count,0.);
        if(dispersion<=1e-9)return std::max(0.,mean-count+(count?count*std::log(count/mean):0));
        double r=1/dispersion;
        return std::max(0.,(count?count*std::log(count/mean):0)-(count+r)*std::log((count+r)/(mean+r)));
    }
    double observe(double count,const Options& options) {
        profile.observe(count);
        if(baseline.size()<options.coverage_window) {baseline.push_back(count);return 0;}
        double alpha=options.coverage_model=="poisson"?0:dispersion();
        // Explicit NB uses a weak overdispersion floor when local MAD is zero.
        if(options.coverage_model=="negative-binomial")alpha=std::max(alpha,.05);
        double value=loss(count,mean(),alpha);penalty+=value;++scored;return value;
    }
    double score() const {return scored?-penalty/scored:0;}
};
struct Stability {
    double within=0,jumps=0,between=0,support=0;
    double score() const {return -(within+jumps+between+support);}
};
inline Stability stability(const Coverage& left,const Coverage& right) {
    const auto& l=left.profile;const auto& r=right.profile;
    Stability s;
    s.within=(l.variance()+r.variance())/2;
    s.jumps=(l.roughness()+r.roughness())/2;
    if(l.samples && r.samples)s.between=(l.mean_log-r.mean_log)*(l.mean_log-r.mean_log);
    // A short/absent flank must not earn a perfect score from zero observations.
    s.support=1./(l.samples+1)+1./(r.samples+1);return s;
}
struct Side {
    Path nodes; // starts at seed endpoint, grows outwards
    Coverage coverage;
    size_t branches=0,since_sample=0;
    std::string stop;
    std::string context="unexamined",context_prefix;
    size_t context_clear_at=0;
};
struct WalkResult {std::vector<Side> paths;size_t expanded=0;bool beam_pruned=false,work_limited=false;};
inline bool better_side(const Side& a,const Side& b) {
    if((a.nodes.size()>1)!=(b.nodes.size()>1))return a.nodes.size()>1;
    if(a.coverage.profile.loss()!=b.coverage.profile.loss())return a.coverage.profile.loss()<b.coverage.profile.loss();
    if(a.nodes.size()!=b.nodes.size())return a.nodes.size()>b.nodes.size();
    return a.nodes<b.nodes;
}

template<class Graph>
WalkResult extend(Graph& graph,const Seed& seed,bool backward,size_t budget,
                  const std::set<Node>& array_edges,const Options& options,const std::set<Node>* shared_context=nullptr,
                  const std::vector<FlankContext>* contexts=nullptr) {
    WalkResult result;
    const auto& first=backward?seed.left:seed.right;
    Node endpoint=backward?seed.path.front():seed.path.back();
    Side initial;initial.nodes={endpoint};
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
        for(const auto& old:result.paths)if(old.nodes==state.nodes && old.stop==state.stop && old.context==state.context)return;
        result.paths.push_back(std::move(state));
        // Bound terminal storage as well as the frontier.
        if(result.paths.size()>options.beam_width*2) {
            std::sort(result.paths.begin(),result.paths.end(),better_side);
            result.paths.resize(options.beam_width);result.beam_pruned=true;
        }
    };
    while(!frontier.empty()) {
        std::sort(frontier.begin(),frontier.end(),better_side);
        Side current=std::move(frontier.front());frontier.erase(frontier.begin());
        std::unordered_set<Node> visited(seed.path.begin(),seed.path.end());
        visited.insert(current.nodes.begin(),current.nodes.end());
        // Nonbranching stretches are walked in place: no full path copy per k-mer.
        while(true) {
            if(current.nodes.size()-1>=budget) {finish(std::move(current),"length_limit");break;}
            if(result.expanded>=options.max_states) {
                result.work_limited=true;finish(std::move(current),"work_limit");
                for(auto& p:frontier)finish(std::move(p),"work_limit");
                frontier.clear();break;
            }
            auto arcs=current.nodes.size()==1?first:graph.neighbors(current.nodes.back(),backward);
            std::vector<Arc> usable;bool cycle=false,masked=false;
            for(auto a:arcs) {
                if(visited.count(a.to)){cycle=true;continue;}
                if(array_edges.count(a.edge)){masked=true;continue;}
                usable.push_back(a);
            }
            std::sort(usable.begin(),usable.end(),[](auto a,auto b){return a.to<b.to;});
            if(usable.empty()){finish(std::move(current),cycle?"cycle":masked?"array_context":"dead_end");break;}
            bool branch=usable.size()>1;
            if(branch)++current.branches;
            std::vector<Side> children;
            for(const auto& a:usable) {
                if(result.expanded>=options.max_states){result.work_limited=true;break;}
                ++result.expanded;
                Side next=branch?current:std::move(current);
                if(next.context=="pending") {
                    auto label=graph.label(a.to);
                    next.context_prefix+=backward?label.front():label.back();
                    auto state=context_state(next.context_prefix,*contexts,options.array_repeat_identity);
                    if(state==ContextState::array) {
                        next.context="array";finish(std::move(next),"source_array_context");continue;
                    }
                    if(state==ContextState::external) {
                        next.context="external";next.context_clear_at=next.nodes.size();
                    }
                }
                bool sample=(!shared_context || !shared_context->count(a.edge)) &&
                    next.context!="pending" &&
                    (next.coverage.baseline.empty() || next.since_sample>=graph.k()+1);
                auto before=next.coverage;
                double penalty=sample?next.coverage.observe(a.count,options):0;
                if(sample)next.since_sample=0;else ++next.since_sample;
                if(penalty>options.max_coverage_penalty) {next.coverage=std::move(before);finish(std::move(next),"coverage_break");continue;}
                next.nodes.push_back(a.to);children.push_back(std::move(next));
            }
            if(!branch && children.size()==1) {
                current=std::move(children.front());visited.insert(current.nodes.back());continue;
            }
            for(auto& c:children)frontier.push_back(std::move(c));
            if(frontier.size()>options.beam_width) {
                std::sort(frontier.begin(),frontier.end(),better_side);
                frontier.resize(options.beam_width);result.beam_pruned=true;
            }
            break;
        }
    }
    std::sort(result.paths.begin(),result.paths.end(),better_side);
    // The same partial path can stop before several incompatible successors.
    std::set<std::pair<Path,std::string>> seen;std::vector<Side> unique;
    for(auto& side:result.paths)if(seen.emplace(side.nodes,side.context).second)unique.push_back(std::move(side));
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
};
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
                                   const std::vector<FlankContext>* left_context=nullptr,const std::vector<FlankContext>* right_context=nullptr) {
    size_t seed_bp=seed.path.size()+graph.k()-1;
    size_t budget=options.max_bp>seed_bp?options.max_bp-seed_bp:0;
    // Each side can use the whole budget; combination trims symmetrically and
    // transfers spare capacity, so a blocked side never erases the other flank.
    auto left=extend(graph,seed,true,budget,array_edges,options,shared_context,left_context);
    auto right=extend(graph,seed,false,budget,array_edges,options,shared_context,right_context);
    std::vector<Candidate> results;std::set<std::tuple<Path,std::string,std::string>> seen;size_t combinations=0;
    bool combination_limit=left.paths.size()*right.paths.size()>options.beam_width*4;
    for(size_t diagonal=0;diagonal<left.paths.size()+right.paths.size();++diagonal)
      for(size_t li=0;li<left.paths.size() && li<=diagonal;++li) {
        size_t ri=diagonal-li;if(ri>=right.paths.size())continue;
        if(combinations++>=options.beam_width*4)continue;
        const auto& l=left.paths[li];const auto& r=right.paths[ri];
        size_t ln=std::min(l.nodes.size()-1,budget/2);
        size_t rn=std::min(r.nodes.size()-1,budget-ln);
        ln=std::min(l.nodes.size()-1,budget-rn);
        Candidate c;c.left=l;c.right=r;
        if(ln+1<l.nodes.size())c.left.stop="length_limit";
        if(rn+1<r.nodes.size())c.right.stop="length_limit";
        c.left_length=ln;c.right_length=rn;
        // Recompute evidence only on emitted flanks, excluding trimmed suffixes.
        auto evidence=[&](const Path& nodes,size_t length,bool backwards,size_t clear_at) {
            Coverage cov;size_t since=0;
            for(size_t i=1;i<=length;++i) {
                if(i<clear_at)continue;
                Node edge=missing;for(auto a:graph.neighbors(nodes[i-1],backwards))if(a.to==nodes[i]){edge=a.edge;break;}
                if(shared_context && shared_context->count(edge)){++since;continue;}
                if(!cov.baseline.empty() && since<graph.k()+1){++since;continue;}
                double count=backwards?graph.count(nodes[i],nodes[i-1]):graph.count(nodes[i-1],nodes[i]);
                cov.observe(count,options);since=0;
            }
            return cov;
        };
        c.left.coverage=evidence(l.nodes,ln,true,l.context_clear_at);c.right.coverage=evidence(r.nodes,rn,false,r.context_clear_at);
        Path{}.swap(c.left.nodes);Path{}.swap(c.right.nodes);
        for(size_t i=ln;i>0;--i)c.path.push_back(l.nodes[i]);
        c.seed_begin=ln;c.seed_end=ln+seed_bp;
        c.path.insert(c.path.end(),seed.path.begin(),seed.path.end());
        for(size_t i=1;i<=rn;++i)c.path.push_back(r.nodes[i]);
        // Repeated k-mers already present in the supplied spacer are allowed;
        // new flank vertices must still be disjoint from seed and opposite flank.
        if(std::set<Node>(c.path.begin(),c.path.end()).size()!=std::set<Node>(seed.path.begin(),seed.path.end()).size()+ln+rn)continue;
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
        c.expanded=left.expanded+right.expanded;
        results.push_back(std::move(c));
    }
    return results;
}
} // namespace mcaat::proto
#endif
