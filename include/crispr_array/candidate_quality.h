#ifndef MCAAT_CANDIDATE_QUALITY_H
#define MCAAT_CANDIDATE_QUALITY_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace mcaat::arrays {

// Global Levenshtein distance, with unit cost for each edit. The common DNA
// lengths fit a single Myers bit vector; longer configured lengths use DP.
inline size_t edit_distance(const std::string& a,const std::string& b) {
    if (a.empty()) return b.size();
    if (b.empty()) return a.size();
    if (a.size()>b.size()) return edit_distance(b,a);
    if (a.size()<=64) {
        std::array<uint64_t,256> equality{};
        for (size_t i=0;i<a.size();++i) equality[static_cast<unsigned char>(a[i])]|=uint64_t(1)<<i;
        uint64_t positive=~uint64_t(0),negative=0,high=uint64_t(1)<<(a.size()-1);
        size_t score=a.size();
        for (unsigned char c:b) {
            uint64_t eq=equality[c],vertical=eq|negative;
            uint64_t horizontal=(((eq&positive)+positive)^positive)|eq;
            uint64_t plus=negative|~(horizontal|positive),minus=positive&horizontal;
            if (plus&high) ++score;
            if (minus&high) --score;
            plus=(plus<<1)|1;minus<<=1;
            positive=minus|~(vertical|plus);negative=plus&vertical;
        }
        return score;
    }
    std::vector<size_t> previous(a.size()+1),current(a.size()+1);
    std::iota(previous.begin(),previous.end(),0);
    for (size_t j=0;j<b.size();++j) {
        current[0]=j+1;
        for (size_t i=0;i<a.size();++i)
            current[i+1]=std::min({current[i]+1,previous[i+1]+1,previous[i]+size_t(a[i]!=b[j])});
        previous.swap(current);
    }
    return previous.back();
}

inline double sequence_similarity(const std::string& a,const std::string& b) {
    const auto length=std::max(a.size(),b.size());
    return length ? 1.-double(edit_distance(a,b))/length : 1.;
}

struct QualityOptions {
    double variant_similarity=.90;
    double max_spacer_similarity=.60;
    double max_repeat_spacer_similarity=.60;
    size_t max_cycles=500; // Exclusive, matching cycle enumeration.
    void validate() const {
        for (auto value:{variant_similarity,max_spacer_similarity,max_repeat_spacer_similarity})
            if (!std::isfinite(value) || value<=0 || value>1)
                throw std::invalid_argument("array similarity limits must be finite fractions in (0, 1]");
        if (max_cycles<3) throw std::invalid_argument("max-cycles must be at least 3");
    }
};

struct SpacerFamily {
    std::string representative;              // An observed sequence, never consensus.
    std::vector<std::string> variants;        // Every distinct observed sequence retained.
};

struct Quality {
    std::vector<SpacerFamily> families;
    size_t variants=0;
    double spacer_similarity=0,repeat_spacer_similarity=0;
    bool similarities_computed=false;
    std::string reason; // Empty means the candidate passed these specified checks.
};

inline Quality assess_quality(const std::string& repeat,std::vector<std::string> spacers,
                              size_t cycle_count,const QualityOptions& options) {
    Quality result;
    std::sort(spacers.begin(),spacers.end());
    spacers.erase(std::unique(spacers.begin(),spacers.end()),spacers.end());
    result.variants=spacers.size();
    if (cycle_count>=options.max_cycles) {result.reason="cycle_limit";return result;}
    if (spacers.size()<2) {result.reason="fewer_than_two_spacers";return result;}
    std::vector<size_t> lengths;for (const auto& s:spacers) lengths.push_back(s.size());
    std::sort(lengths.begin(),lengths.end());
    double median=lengths.size()%2?lengths[lengths.size()/2]:
        (double(lengths[lengths.size()/2-1])+lengths[lengths.size()/2])/2;
    if (median<.5*repeat.size()) {result.reason="spacer_repeat_length_ratio";return result;}
    const size_t n=spacers.size();
    std::vector<size_t> distances(n*n);
    for (size_t i=0;i<n;++i) for (size_t j=0;j<i;++j)
        distances[i*n+j]=distances[j*n+i]=edit_distance(spacers[i],spacers[j]);
    // Deterministic complete-link groups: every pair in a group meets the
    // threshold. A-B and B-C similarity alone must not merge dissimilar A-C.
    std::vector<std::vector<size_t>> groups;
    for (size_t i=0;i<n;++i) {
        size_t chosen=groups.size();
        for (size_t c=0;c<groups.size();++c) {
            bool fits=true;
            for (size_t j:groups[c])
                if (double(distances[i*n+j])> (1.-options.variant_similarity)*std::max(spacers[i].size(),spacers[j].size())+1e-9) {fits=false;break;}
            if (fits) {chosen=c;break;}
        }
        if (chosen==groups.size()) groups.push_back({});
        groups[chosen].push_back(i);
    }
    for (const auto& group:groups) {
        size_t best=group.front(),best_sum=std::numeric_limits<size_t>::max();
        for (size_t i:group) {
            size_t sum=0;for (size_t j:group) sum+=distances[i*n+j];
            if (sum<best_sum) {best=i;best_sum=sum;}
        }
        SpacerFamily family;family.representative=spacers[best];
        for (size_t i:group) family.variants.push_back(spacers[i]);
        result.families.push_back(std::move(family));
    }
    std::sort(result.families.begin(),result.families.end(),[](const auto& a,const auto& b) {return a.representative<b.representative;});
    if (result.families.size()<2) {result.reason="fewer_than_two_spacer_families";return result;}
    size_t pairs=0;
    for (size_t i=0;i<result.families.size();++i) {
        const auto& sequence=result.families[i].representative;
        result.repeat_spacer_similarity+=sequence_similarity(repeat,sequence);
        for (size_t j=0;j<i;++j) {
            result.spacer_similarity+=sequence_similarity(sequence,result.families[j].representative);++pairs;
        }
    }
    result.spacer_similarity/=pairs;
    result.repeat_spacer_similarity/=result.families.size();
    result.similarities_computed=true;
    if (result.spacer_similarity>options.max_spacer_similarity+1e-12)
        result.reason="similar_spacer_families";
    else if (result.repeat_spacer_similarity>options.max_repeat_spacer_similarity+1e-12)
        result.reason="repeat_like_spacers";
    return result;
}
} // namespace mcaat::arrays
#endif
