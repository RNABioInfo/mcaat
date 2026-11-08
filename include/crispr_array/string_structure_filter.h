#ifndef MCAAT_STRING_STRUCTURE_FILTER_H
#define MCAAT_STRING_STRUCTURE_FILTER_H
#include "crispr_array/candidate_quality.h"
#include "spoa/spoa.hpp"
#include <map>

namespace mcaat::arrays {
struct StringFilterOptions {
    bool enabled=true;
    double column_support=.85, motif_identity=.90, motif_support=.85;
    double conserved_spacer_fraction=.80, max_contrast=.15;
    void validate() const {
        for (double value:{column_support,motif_identity,motif_support,conserved_spacer_fraction})
            if (!std::isfinite(value) || value<=0 || value>1)
                throw std::invalid_argument("string-filter fractions must be finite and in (0, 1]");
        if (!std::isfinite(max_contrast) || max_contrast<0 || max_contrast>1)
            throw std::invalid_argument("string-max-contrast must be finite and in [0, 1]");
    }
};
inline bool string_filter_option(std::string key) {
    if (key.rfind("--",0)==0)key.erase(0,2);
    return key=="string-filter" || key=="string-column-support" || key=="string-motif-identity" ||
           key=="string-motif-support" || key=="string-conserved-spacer-fraction" || key=="string-max-contrast";
}
inline bool set_string_filter_option(std::string key,const std::string& value,StringFilterOptions& options) {
    if (key.rfind("--",0)==0)key.erase(0,2);
    if (!string_filter_option(key))return false;
    if (key=="string-filter") {
        if (value!="on" && value!="off")throw std::invalid_argument("string-filter must be on or off");
        options.enabled=value=="on";return true;
    }
    size_t used=0;double number=0;
    try {number=std::stod(value,&used);}catch (...) {used=0;}
    if (used!=value.size() || !used || !std::isfinite(number) || number>1 ||
        (key=="string-max-contrast"?number<0:number<=0))
        throw std::invalid_argument(key+" has an invalid fraction");
    if (key=="string-column-support")options.column_support=number;
    else if (key=="string-motif-identity")options.motif_identity=number;
    else if (key=="string-motif-support")options.motif_support=number;
    else if (key=="string-conserved-spacer-fraction")options.conserved_spacer_fraction=number;
    else options.max_contrast=number;
    return true;
}
struct StringEvidence {
    bool evaluated=false, alignment_scored=false, alignment_skipped=false;
    bool short_motif=false, conserved_spacer=false;
    size_t observations=0,families=0,repeat_columns=0,spacer_columns=0;
    double motif_fraction=0,mean_motif_identity=0;
    double repeat_conservation=0,spacer_conservation=0,conserved_spacer_fraction=0;
    std::string reason() const {
        if (short_motif && conserved_spacer)return "short_motif_across_unit+conservation_across_spacer";
        if (short_motif)return "short_motif_across_unit";
        if (conserved_spacer)return "conservation_across_spacer";
        return {};
    }
};
inline size_t dna_base(char c) {
    switch(c) {case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;}
    throw std::invalid_argument("string filter requires A/C/G/T observations");
}
// Best substitution-only periodic template, with at least three nominal copies.
// No template or consensus replaces any input sequence.
inline double short_motif_identity(const std::string& sequence) {
    std::vector<size_t> bases;for(char c:sequence)bases.push_back(dna_base(c));
    size_t best=0;
    for(size_t period=1;period<=bases.size()/3;++period) {
        std::vector<std::array<size_t,4>> counts(period);
        for(size_t i=0;i<bases.size();++i)++counts[i%period][bases[i]];
        size_t matches=0;for(const auto& column:counts)matches+=*std::max_element(column.begin(),column.end());
        best=std::max(best,matches);
    }
    return bases.empty()?0:double(best)/bases.size();
}
inline std::vector<std::string> observed_alignment(const std::vector<std::string>& sequences) {
    if(sequences.empty())throw std::invalid_argument("cannot align empty observation set");
    std::map<std::string,size_t> unique;
    for(const auto& sequence:sequences) {
        if(sequence.empty())throw std::invalid_argument("cannot align empty sequence");
        for(char c:sequence)dna_base(c);
        unique.emplace(sequence,0);
    }
    if(unique.size()==1)return sequences;
    spoa::Graph graph;
    auto engine=spoa::AlignmentEngine::Create(spoa::AlignmentType::kNW,3,-5,-3);
    size_t index=0;
    // Canonical insertion order makes alignment independent of graph IDs and
    // input row order. Duplicate strings do not alter the alignment topology.
    for(auto& [sequence,id]:unique) {
        id=index++;
        graph.AddAlignment(engine->Align(sequence,graph),sequence);
    }
    auto aligned=graph.GenerateMultipleSequenceAlignment(false);
    if(aligned.size()!=unique.size())throw std::runtime_error("alignment changed observation count");
    for(const auto& [sequence,id]:unique) {
        auto ungapped=aligned.at(id);ungapped.erase(std::remove(ungapped.begin(),ungapped.end(),'-'),ungapped.end());
        if(ungapped!=sequence || aligned[id].size()!=aligned.front().size())
            throw std::runtime_error("alignment changed an observed string");
    }
    std::vector<std::string> result;result.reserve(sequences.size());
    for(const auto& sequence:sequences)result.push_back(aligned.at(unique.at(sequence)));
    return result;
}
inline std::vector<double> base_conservation(const std::vector<std::string>& aligned,
                                              const std::vector<double>& weights) {
    if(aligned.empty() || aligned.size()!=weights.size())throw std::invalid_argument("invalid alignment weights");
    double total=0;for(double w:weights) {
        if(!std::isfinite(w) || w<=0)throw std::invalid_argument("invalid alignment weight");
        total+=w;
    }
    std::vector<std::array<double,4>> counts(aligned.front().size());
    for(size_t row=0;row<aligned.size();++row) {
        if(aligned[row].size()!=counts.size())throw std::invalid_argument("unequal alignment widths");
        for(size_t col=0;col<counts.size();++col)
            if(aligned[row][col]!='-')counts[col][dna_base(aligned[row][col])]+=weights[row]/total;
    }
    std::vector<double> result;
    for(const auto& column:counts)result.push_back(*std::max_element(column.begin(),column.end()));
    return result;
}
inline StringEvidence assess_strings(const std::vector<std::pair<std::string,std::string>>& observations,
                                      const Quality& quality,const StringFilterOptions& options,
                                      size_t alignment_limit) {
    StringEvidence result;
    if(!options.enabled)return result;
    options.validate();
    if(observations.empty() || quality.families.empty())throw std::invalid_argument("missing string-filter observations/families");
    result.evaluated=true;result.observations=observations.size();result.families=quality.families.size();
    std::map<std::string,size_t> family;
    for(size_t i=0;i<quality.families.size();++i)
        for(const auto& spacer:quality.families[i].variants)
            if(!family.emplace(spacer,i).second)throw std::runtime_error("spacer belongs to multiple families");
    std::vector<size_t> counts(quality.families.size());
    for(const auto& row:observations)++counts.at(family.at(row.second));
    std::vector<double> weights;std::vector<std::string> repeats,spacers;
    for(const auto& [repeat,spacer]:observations) {
        const double weight=1./result.families/counts.at(family.at(spacer));
        weights.push_back(weight);repeats.push_back(repeat);spacers.push_back(spacer);
        const double identity=short_motif_identity(repeat+spacer);
        result.mean_motif_identity+=weight*identity;
        if(identity+1e-12>=options.motif_identity)result.motif_fraction+=weight;
    }
    result.short_motif=result.families>=2 && result.motif_fraction+1e-12>=options.motif_support;
    // Keep the existing bounded-work policy for merged candidates at the cycle
    // cap. Report skipped alignment explicitly; motif scoring is still performed.
    if(observations.size()>=alignment_limit) {result.alignment_skipped=true;return result;}
    auto rp=base_conservation(observed_alignment(repeats),weights);
    auto sp=base_conservation(observed_alignment(spacers),weights);
    result.alignment_scored=true;result.repeat_columns=rp.size();result.spacer_columns=sp.size();
    result.repeat_conservation=std::accumulate(rp.begin(),rp.end(),0.)/rp.size();
    result.spacer_conservation=std::accumulate(sp.begin(),sp.end(),0.)/sp.size();
    result.conserved_spacer_fraction=double(std::count_if(sp.begin(),sp.end(),[&](double x) {
        return x+1e-12>=options.column_support;
    }))/sp.size();
    result.conserved_spacer=result.families>=3 &&
        result.conserved_spacer_fraction+1e-12>=options.conserved_spacer_fraction &&
        result.repeat_conservation-result.spacer_conservation<=options.max_contrast+1e-12;
    return result;
}
} // namespace mcaat::arrays
#endif
