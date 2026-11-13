#ifndef MCAAT_PROTO_PATH_QUALITY_H
#define MCAAT_PROTO_PATH_QUALITY_H
#include "protospacer_detection/extension_engine.h"
#include <array>

namespace mcaat::proto {
// Reference-free evidence checks. Scores are heuristics, not probabilities of
// correctness. Adjacent k-mer counts are correlated; no independence is assumed.
struct FlankQuality {
    size_t bases=0,observations=0,windows=0,low_windows=0,revisited=0,kmers=0;
    size_t entrance_observations=0,downstream_observations=0,transition_pairs=0;
    double singleton_fraction=0,median_count=0,block_variation=0,largest_step=0,loss=0;
    double entrance_energy=0,downstream_energy=0,local_roughness=0,branch_loss=0;
    bool depth_break=false;
    std::vector<std::string> reasons;
};
inline FlankQuality flank_quality(const std::string& sequence,size_t k,const std::vector<double>& counts,
        const std::vector<size_t>& positions={},const std::vector<bool>& entrance={},
        const std::vector<double>& competing_max={}) {
    if((!positions.empty()&&positions.size()!=counts.size()) ||
       (!entrance.empty()&&entrance.size()!=counts.size()) ||
       (!competing_max.empty()&&competing_max.size()!=counts.size()))
        throw std::runtime_error("inconsistent path-quality observations");
    FlankQuality q;q.bases=sequence.size();q.observations=counts.size();
    for(size_t i=0;i+64<=sequence.size();i+=16) {
        ++q.windows;q.low_windows+=triplet_entropy(sequence,i,64)<3.;
    }
    std::set<std::string> words;
    for(size_t i=0;i+k<=sequence.size();++i){++q.kmers;q.revisited+=!words.insert(sequence.substr(i,k)).second;}
    q.median_count=median(counts);
    for(auto c:counts)q.singleton_fraction+=c<=1;
    if(!counts.empty())q.singleton_fraction/=counts.size();
    // Keep the former block summaries as diagnostics. Global flatness does not
    // distinguish a smooth abundance trend from a wrong, flat continuation.
    std::vector<double> blocks;
    for(size_t i=0;i+64<=counts.size();i+=64)
        blocks.push_back(median(std::vector<double>(counts.begin()+i,counts.begin()+i+64)));
    double center=median(blocks);
    for(size_t i=0;i<blocks.size();++i) {
        double d=std::log2((1+blocks[i])/(1+center));q.block_variation+=d*d;
        if(i) {
            double lo=std::max(1.,std::min(blocks[i-1],blocks[i])),hi=std::max(blocks[i-1],blocks[i]);
            q.largest_step=std::max(q.largest_step,std::log2((1+hi)/(1+lo)));
            if(hi>=4*lo && hi-lo>=3*std::sqrt(hi+lo))q.depth_break=true;
        }
    }
    if(!blocks.empty())q.block_variation/=blocks.size();
    // The expected source-overlap interval marks the entrance, but no longer
    // deletes observations outside the actual source path. Do not fabricate
    // an adjacent transition across omitted source edges.
    for(size_t i=0;i<counts.size();++i) {
        bool entering=!entrance.empty()&&entrance[i];
        if(entering)++q.entrance_observations;
        else {
            ++q.downstream_observations;
            if(!competing_max.empty())q.branch_loss+=std::log2(
                (1+std::max(counts[i],competing_max[i]))/(1+counts[i]));
        }
        if(!i || (!positions.empty()&&positions[i]!=positions[i-1]+1))continue;
        double delta=std::log2((1+counts[i])/(1+counts[i-1]));
        ++q.transition_pairs;
        if(entering || (!entrance.empty()&&entrance[i-1]))q.entrance_energy+=delta*delta;
        else q.downstream_energy+=delta*delta;
    }
    if(q.transition_pairs)q.local_roughness=(q.entrance_energy+q.downstream_energy)/q.transition_pairs;
    if(q.downstream_observations)q.branch_loss/=q.downstream_observations;
    if(q.windows && double(q.low_windows)/q.windows>.2)q.reasons.push_back("low_complexity");
    if(q.kmers && double(q.revisited)/q.kmers>.2)q.reasons.push_back("repeat_dominated");
    if(q.downstream_observations<64)q.reasons.push_back("insufficient_downstream_evidence");
    if(q.singleton_fraction>.5)q.reasons.push_back("singleton_dominated");
    // A jump is evidence to rank, not a hard rejection: valid source exits can
    // contain large changes and low-count bridges. Exact source counts never
    // contribute to this score or its baseline.
    // Do not penalize an absent second flank as though it were observed bad
    // sequence. Score each existing flank, then use the worse flank below.
    q.loss=q.local_roughness+q.branch_loss+q.singleton_fraction+
        (q.windows?double(q.low_windows)/q.windows:0)+(q.kmers?double(q.revisited)/q.kmers:0)+1./(1+counts.size());
    return q;
}
struct PathQuality {
    FlankQuality left,right;
    double score=0;
    std::vector<std::string> reasons;
    bool accepted() const{return reasons.empty();}
    std::string reason_text()const {
        std::string s;for(const auto& r:reasons){if(!s.empty())s+=',';s+=r;}return s.empty()?"pass":s;
    }
};
// Scores now measure local changes, on a different scale from block variance.
// A margin of .1 means at most 10% worse loss, not an absolute .1 score gap.
inline bool near_quality_score(double best,double candidate,double relative_margin) {
    return candidate>=best-relative_margin*std::max(1e-12,-best)-1e-12;
}
// Compare branches on the same flank. Opposite flanks and shorter prefixes do
// not by themselves conflict. Shared (k+1)-mers allow small bubbles/indels to
// count as the same route; this is not a nucleotide-identity calculation.
inline bool distinct_extension_routes(const Candidate& a,const Candidate& b,size_t k) {
    for(bool back:{true,false}) {
        auto x=back?a.sequence.substr(0,a.seed_begin):a.sequence.substr(a.seed_end);
        auto y=back?b.sequence.substr(0,b.seed_begin):b.sequence.substr(b.seed_end);
        if(std::min(x.size(),y.size())<std::max(size_t(128),k+1))continue;
        std::set<std::string> one,two;
        for(size_t i=0;i+k<x.size();++i)one.insert(x.substr(i,k+1));
        for(size_t i=0;i+k<y.size();++i)two.insert(y.substr(i,k+1));
        size_t common=0;for(const auto& word:one)common+=two.count(word);
        if(double(common)<.8*std::min(one.size(),two.size()))return true;
    }
    return false;
}
template<class Graph>
PathQuality assess_path_quality(Graph& graph,const Candidate& c,const std::set<Node>& shared,
                               SourceCycleGuard& source_guard,size_t minimum_flank) {
    PathQuality q;
    if(std::max(c.seed_begin,c.sequence.size()-c.seed_end)<minimum_flank)q.reasons.push_back("requested_flank_not_reached");
    for(bool back:{true,false}) {
        auto sequence=back?c.sequence.substr(0,c.seed_begin):c.sequence.substr(c.seed_end);
        if(sequence.empty())continue;
        auto& side=back?c.left:c.right;
        Path path=back?Path(c.path.begin(),c.path.begin()+c.seed_begin+1):Path(c.path.begin()+c.seed_end-graph.k(),c.path.end());
        if(back){std::reverse(sequence.begin(),sequence.end());std::reverse(path.begin(),path.end());}
        std::vector<double> counts,competing;
        std::vector<size_t> positions;
        std::vector<bool> entrance;
        const size_t entrance_end=std::max({graph.k()+1,side.context_clear_at,side.source_mask_until});
        for(size_t i=1;i<path.size();++i) {
            if(i<graph.k()+1)continue;
            const auto choices=graph.neighbors(path[i-1],back);
            double maximum=0;for(const auto& arc:choices)maximum=std::max(maximum,arc.count);
            for(const auto& arc:choices)if(arc.to==path[i]) {
                if(!shared.count(arc.edge) && !source_guard.contains(graph,path[i-1],arc,back)) {
                    counts.push_back(arc.count);positions.push_back(i);
                    entrance.push_back(i<entrance_end);competing.push_back(maximum);
                }
                break;
            }
        }
        auto& out=back?q.left:q.right;out=flank_quality(sequence,graph.k(),counts,positions,entrance,competing);
        for(const auto& reason:out.reasons)q.reasons.push_back(std::string(back?"left:":"right:")+reason);
    }
    q.score=-std::max(q.left.loss,q.right.loss);return q;
}
} // namespace mcaat::proto
#endif
