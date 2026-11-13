#ifndef MCAAT_PROTO_RETAINED_TARGET_H
#define MCAAT_PROTO_RETAINED_TARGET_H
#include "protospacer_detection/extension_engine.h"
#include "protospacer_detection/target_placement.h"
#include "protospacer_detection/path_quality.h"

namespace mcaat::proto {
// All strings are oriented like the supplied, immutable protospacer.
struct TargetContext {
    std::string spacer,repeat;
    size_t begin=0,end=0;
};
struct PlacementCheck {
    std::string core,phase,reason,sequence;
    double score=0;
    size_t query_offset=0,query_begin=0,query_end=0,expanded=0;
    bool extended=false,limited=false;
};
struct QualityCheck {
    std::string sequence;
    size_t seed_begin=0,seed_end=0;
    PathQuality quality;
};
struct TargetExtension {
    std::string status="not_extended",reason;
    Candidate selected;
    std::vector<Candidate> alternatives; // one placement, at most twice the beam width after flank splitting
    bool has_extension=false,ambiguous=false,search_limited=false;
    size_t compared=0,expanded=0;
    double runner_score=0;
    bool has_runner=false;
    SearchDiagnostics diagnostics;
    std::set<std::string> stops;
    TargetPlacement placement;
    size_t placement_work=0,placement_candidates=0;
    bool variant_placement=false,placement_limited=false;
    std::vector<PlacementCheck> placement_checks;
    std::vector<QualityCheck> quality_checks;
    bool quality_unresolved=false;
    size_t competing_quality_routes=0;
};
// Compare against one fixed best-score threshold. A pairwise "near tie"
// comparator would be non-transitive and cannot safely be passed to std::sort.
// Within the existing ambiguity margin, prefer more validated sequence rather
// than rewarding a prematurely stopped, slightly flatter branch.
inline double select_supported_near_tie(std::vector<Candidate>& candidates,const Options& options) {
    const double best=candidates.front().score;
    const auto original=stability(candidates.front().left.coverage,candidates.front().right.coverage);
    const double ambiguity_evidence=best+original.support;
    const double within_limit=original.within+original.jumps+original.weak+options.prediction_margin;
    if(options.resolve_graph_ties || options.prediction_margin<=0)return ambiguity_evidence;
    size_t selected=candidates.size();
    for(size_t i=0;i<candidates.size();++i) {
        // Keep the support term in selection: removing it can make an absent
        // flank look better and even shorten the original winning prediction.
        if(candidates[i].score<best-options.prediction_margin-1e-12)continue;
        const auto s=stability(candidates[i].left.coverage,candidates[i].right.coverage);
        // Better agreement of the two mean depths must not compensate for
        // substantially less stable sequence within an individual flank.
        if(s.within+s.jumps+s.weak>within_limit+1e-12)continue;
        if(selected==candidates.size() || candidates[i].sequence.size()>candidates[selected].sequence.size())selected=i;
    }
    if(selected && selected<candidates.size())std::rotate(candidates.begin(),candidates.begin()+selected,candidates.begin()+selected+1);
    return ambiguity_evidence;
}
template<class Graph>
void refresh_target(Graph& graph,Candidate& c,const std::vector<FlankContext>& left_models,
                    const std::vector<FlankContext>& right_models,const std::set<Node>& shared,const Options& options) {
    auto left=c.sequence.substr(0,c.seed_begin);std::reverse(left.begin(),left.end());
    auto right=c.sequence.substr(c.seed_end);
    size_t minimum=options.min_external_flank?options.min_external_flank:graph.k()+1;
    auto check_side=[&](const std::string& seq,const auto& models,bool backward) {
        if(seq.empty())return;
        auto state=context_state(seq,models,options.array_repeat_identity);
        if((state!=ContextState::external && !(options.allow_source_array_crossing && state==ContextState::array)) || seq.size()<minimum) {
            if(backward)trim_left(c);else trim_right(c,graph.k());
            auto& side=backward?c.left:c.right;
            side.context=state==ContextState::array?"array":"incomplete";
            side.stop=state==ContextState::array?"source_array_context":
                state==ContextState::pending?"incomplete_context":"insufficient_external_sequence";
        }
    };
    check_side(left,left_models,true);check_side(right,right_models,false);
    std::vector<FlankContext> models=left_models;models.insert(models.end(),right_models.begin(),right_models.end());
    SourceCycleGuard source_guard(graph.k(),&models);
    c.external_left=c.external_right=0;
    for(size_t i=1;i<c.path.size();++i) {
        if(shared.count(graph.edge(c.path[i-1],c.path[i])) ||
           (options.allow_source_array_crossing && source_guard.contains(graph,c.path[i-1],c.path[i],false)))continue;
        if(i<=c.seed_begin && c.seed_begin-i+1>=c.left.context_clear_at)++c.external_left;
        if(i>c.seed_end-graph.k() && i-(c.seed_end-graph.k())>=c.right.context_clear_at)++c.external_right;
    }
    if(c.left_length && (!c.external_left || (c.left.crosses_source_array && c.external_left<minimum))) {
        auto reason=c.external_left?"insufficient_external_sequence":"array_context_only";
        trim_left(c);c.left.stop=reason;
    }
    if(c.right_length && (!c.external_right || (c.right.crosses_source_array && c.external_right<minimum))) {
        auto reason=c.external_right?"insufficient_external_sequence":"array_context_only";
        trim_right(c,graph.k());c.right.stop=reason;
    }
    Path l(c.path.begin(),c.path.begin()+c.seed_begin+1);std::reverse(l.begin(),l.end());
    Path r(c.path.begin()+c.seed_end-graph.k(),c.path.end());
    c.left.coverage=path_coverage(graph,l,l.size()-1,true,c.left.context_clear_at,&shared,options,&source_guard);
    c.right.coverage=path_coverage(graph,r,r.size()-1,false,c.right.context_clear_at,&shared,options,&source_guard);
    c.score=stability(c.left.coverage,c.right.coverage).score();
}
template<class Graph>
TargetExtension extend_retained_target(Graph& graph,const std::string& sequence,
                                      const std::vector<TargetContext>& contexts,const Options& options,
                                      const std::set<std::string>* source_cycles=nullptr,bool exact_source_interval=true) {
    TargetExtension result;
    result.placement=score_target_placement(sequence,sequence,0);
    if(sequence.size()<graph.k()){result.reason="shorter_than_graph_k";return result;}
    if(sequence.size()>options.max_bp){result.reason="seed_exceeds_length_cap";return result;}
    if(contexts.empty()){result.reason="source_context_missing";return result;}
    if(contexts.size()>1024){result.reason="source_context_limit";result.search_limited=true;return result;}
    if(AdapterGuard(options.adapters).contains(sequence)){result.reason="adapter_in_seed";return result;}
    auto path=graph.trace(sequence);
    if(path.empty()){result.reason="seed_not_mapped";return result;}
    std::set<Node> shared;
    std::vector<FlankContext> left_models,right_models;
    for(const auto& ctx:contexts) {
        if(ctx.end<ctx.begin || ctx.end>ctx.spacer.size() || ctx.end-ctx.begin!=sequence.size() ||
                (exact_source_interval && ctx.spacer.substr(ctx.begin,ctx.end-ctx.begin)!=sequence))
            throw std::runtime_error("retained target context does not contain its seed");
        auto edges=graph.array_context(ctx.repeat,ctx.spacer);
        shared.insert(edges.first.begin(),edges.first.end());
        auto left=flank_contexts(ctx.spacer,ctx.begin,ctx.end,{ctx.repeat},true);
        auto right=flank_contexts(ctx.spacer,ctx.begin,ctx.end,{ctx.repeat},false);
        left_models.insert(left_models.end(),left.begin(),left.end());
        right_models.insert(right_models.end(),right.begin(),right.end());
    }
    if(source_cycles)for(const auto& cycle:*source_cycles) {
        left_models.push_back({"","",cycle});right_models.push_back({"","",cycle});
    }
    Seed seed{path,0,sequence.size(),graph.neighbors(path.front(),true),graph.neighbors(path.back(),false)};
    auto local=options;local.max_states=std::min(options.max_states,options.max_query_states/2);
    if(!local.max_states){result.reason="query_work_limit";result.search_limited=true;return result;}
    auto paths=reconstruct(graph,seed,{},local,&shared,&left_models,&right_models,&result.diagnostics);
    std::vector<Candidate> candidates;std::set<Path> seen;
    for(auto& c:paths) {
        result.expanded=std::max(result.expanded,c.expanded);result.search_limited|=c.search_limited;
        refresh_target(graph,c,left_models,right_models,shared,options);
        for(auto& refined:refine_flanks(std::move(c),graph.k(),options)) {
            result.stops.insert(refined.left.stop);result.stops.insert(refined.right.stop);
            if(!(refined.external_left || refined.external_right) || !seen.insert(refined.path).second)continue;
            candidates.push_back(std::move(refined));
        }
    }
    result.compared=candidates.size();
    if(candidates.empty()){result.reason=result.search_limited?"search_limit_no_extension":"no_supported_extension";return result;}
    std::sort(candidates.begin(),candidates.end(),[&](const auto& a,const auto& b){
        return better_prediction(a,sequence.size(),b,sequence.size(),options);
    });
    double evidence=select_supported_near_tie(candidates,options);
    const auto& winner=candidates.front();
    size_t keep_left=winner.seed_begin,keep_right=winner.sequence.size()-winner.seed_end;
    for(size_t i=1;i<candidates.size();++i) {
        const auto& c=candidates[i];
        if(!result.has_runner){result.has_runner=true;result.runner_score=c.score;}
        double other=c.score+stability(c.left.coverage,c.right.coverage).support;
        if(other<evidence-options.prediction_margin-1e-12)continue;
        result.ambiguous=true;
        if(options.resolve_graph_ties) {
            auto common=common_flanks(winner,0,c,0);
            keep_left=std::min(keep_left,common.left);keep_right=std::min(keep_right,common.right);
        }
    }
    if(options.report_alternatives)result.alternatives=candidates;
    result.selected=std::move(candidates.front());
    if(options.resolve_graph_ties && (keep_left<result.selected.seed_begin || keep_right<result.selected.sequence.size()-result.selected.seed_end)) {
        retain_flanks(result.selected,keep_left,keep_right,graph.k());result.selected.ambiguity_trimmed=true;
        refresh_target(graph,result.selected,left_models,right_models,shared,options);
    }
    if(!(result.selected.external_left || result.selected.external_right)) {
        result.reason="ambiguous_extension_withheld";return result;
    }
    if(result.selected.sequence.substr(result.selected.seed_begin,sequence.size())!=sequence)
        throw std::runtime_error("extension changed a supplied protospacer");
    result.has_extension=true;result.status=result.ambiguous?"extended_ambiguous":"extended";
    result.reason=result.search_limited?"bounded_search":result.ambiguous?"best_candidate_near_tie":"best_candidate";
    return result;
}

// Keep an existing extension. If exact placement produces none, search nearby
// graph alleles instead of forcing the query sequence into a viral target.
// Original input strings and source-array models are never edited.
template<class Graph>
TargetExtension extend_target_with_fallback(Graph& graph,const std::string& query,
        const std::vector<TargetContext>& contexts,const Options& options,const std::set<std::string>* source_cycles=nullptr,
        const TargetExtension* previous_attempt=nullptr) {
    auto exact=previous_attempt?*previous_attempt:extend_retained_target(graph,query,contexts,options,source_cycles);
    if(exact.has_extension || !options.target_variant_fallback || contexts.empty() || query.size()<15 ||
            contexts.size()>1024 || query.size()>options.max_bp ||
            AdapterGuard(options.adapters).contains(query))return exact;
    size_t extra=graph.k()>query.size()?graph.k()-query.size():0;
    size_t available_left=std::numeric_limits<size_t>::max(),available_right=available_left;
    for(const auto& ctx:contexts) {
        if(ctx.end<ctx.begin || ctx.end>ctx.spacer.size() || ctx.spacer.substr(ctx.begin,ctx.end-ctx.begin)!=query)
            throw std::runtime_error("variant source context does not contain the original input");
        available_left=std::min(available_left,ctx.begin);available_right=std::min(available_right,ctx.spacer.size()-ctx.end);
    }
    size_t left_padding=extra>available_right?extra-available_right:0;
    if(left_padding>available_left){exact.reason="short_query_context_insufficient";return exact;}
    auto placements=previous_attempt?
        place_target_masked(graph,query,options.max_target_placement_states-exact.placement_work,options.max_seeds):
        place_target(graph,query,left_padding,options.max_target_placement_states,options.max_seeds,options.path_quality_filter);
    if(previous_attempt) {
        std::set<std::string> checked{query};for(const auto& check:exact.placement_checks)checked.insert(check.core);
        placements.values.erase(std::remove_if(placements.values.begin(),placements.values.end(),[&](const auto& p){return checked.count(p.sequence);}),placements.values.end());
    }
    exact.placement_work+=placements.work;exact.placement_candidates+=placements.values.size();exact.placement_limited|=placements.limited;
    auto adjusted=contexts;for(auto& ctx:adjusted){ctx.begin-=left_padding;ctx.end+=extra-left_padding;}
    size_t remaining=options.max_query_states>exact.expanded?options.max_query_states-exact.expanded:0;
    std::vector<TargetExtension> alternatives;
    std::vector<PlacementCheck> checks=exact.placement_checks;
    auto evaluate=[&](const TargetPlacement& placement,Options local,const std::string& phase) {
        local.left_seed_padding=placement.query_offset;
        local.right_seed_padding=placement.sequence.size()-placement.query_offset-query.size();
        auto candidate=extend_retained_target(graph,placement.sequence,adjusted,local,source_cycles,false);
        remaining-=std::min(remaining,candidate.expanded);
        exact.placement_limited|=candidate.search_limited;
        candidate.placement=placement;candidate.variant_placement=true;
        if(candidate.has_extension) {
            auto penalty=(1-placement.coverage)+(1-placement.identity);
            candidate.selected.score-=penalty;
            if(candidate.has_runner)candidate.runner_score-=penalty;
        }
        checks.push_back({placement.sequence,phase,candidate.reason,candidate.has_extension?candidate.selected.sequence:"",
                          candidate.has_extension?candidate.selected.score:0,placement.query_offset,placement.query_begin,placement.query_end,
                          candidate.expanded,candidate.has_extension,candidate.search_limited});
        return candidate;
    };
    for(size_t i=0;i<placements.values.size();++i) {
        const auto& placement=placements.values[i];
        if(placement.sequence==query)continue;
        if(remaining<2){exact.placement_limited=true;break;}
        auto local=options;
        if(extra>=local.max_flank_bp)continue;
        // Give every placement an initial look at both external flanks before
        // a repetitive early placement can spend the entire query budget.
        local.max_flank_bp=std::min(size_t(75),local.max_flank_bp);
        local.max_query_states=std::max(size_t(2),std::min(size_t(1000),remaining/(placements.values.size()-i)));
        auto candidate=evaluate(placement,local,previous_attempt?"masked_probe":"probe");
        if(!candidate.has_extension)continue;
        alternatives.push_back(std::move(candidate));
    }
    auto better=[&](const auto& a,const auto& b){
        return better_prediction(a.selected,a.placement.query_end-a.placement.query_begin,
                                 b.selected,b.placement.query_end-b.placement.query_begin,options);
    };
    std::sort(alternatives.begin(),alternatives.end(),better);
    const size_t refine=std::min(options.beam_width,alternatives.size());
    for(size_t i=0;i<refine;++i) {
        if(remaining<2){exact.placement_limited=true;break;}
        auto local=options;
        local.max_query_states=std::max(size_t(2),remaining/(refine-i));
        auto candidate=evaluate(alternatives[i].placement,local,previous_attempt?"masked_extend":"extend");
        if(candidate.has_extension)alternatives[i]=std::move(candidate);
    }
    exact.placement_limited|=alternatives.size()>refine;
    alternatives.resize(refine);
    exact.expanded=options.max_query_states-remaining;
    exact.search_limited|=exact.placement_limited;
    if(alternatives.empty()) {
        exact.placement_checks=std::move(checks);
        exact.reason=exact.placement_limited?"variant_search_limit_no_extension":"no_supported_variant_extension";
        if(!options.path_quality_filter&&!previous_attempt&&query.size()>=graph.k()&&remaining>=2&&exact.placement_work<options.max_target_placement_states) {
            auto rescue=options;rescue.noise_adjusted_stability=true;
            return extend_target_with_fallback(graph,query,contexts,rescue,source_cycles,&exact);
        }
        return exact;
    }
    std::sort(alternatives.begin(),alternatives.end(),better);
    auto result=std::move(alternatives.front());
    if(alternatives.size()>1) {
        if(!result.has_runner || alternatives[1].selected.score>result.runner_score)result.runner_score=alternatives[1].selected.score;
        result.has_runner=true;
        if(alternatives[1].selected.score>=result.selected.score-options.prediction_margin) {
            result.ambiguous=true;
            if(options.resolve_graph_ties) {
                exact.reason="ambiguous_variant_placement_withheld";exact.ambiguous=true;
                exact.placement_checks=std::move(checks);return exact;
            }
        }
    }
    result.placement_work=exact.placement_work;result.placement_candidates=exact.placement_candidates;result.placement_limited=exact.placement_limited;
    result.expanded=options.max_query_states-remaining;
    result.placement_checks=std::move(checks);
    result.search_limited|=result.placement_limited;
    result.status=result.ambiguous?"extended_variant_ambiguous":"extended_variant";
    result.reason="graph_variant_placement";return result;
}
// Apply after placement/search, never during short placement probes. The input
// sequence stays in Protospacers.fasta even when every extension is withheld.
template<class Graph>
void filter_target_extension(Graph& graph,TargetExtension& result,const std::string& query,
        const std::vector<TargetContext>& contexts,const std::set<std::string>& cycles,const Options& options) {
    if(!options.path_quality_filter || !result.has_extension)return;
    std::set<Node> shared;
    std::vector<FlankContext> source_models;
    for(const auto& ctx:contexts) {
        auto edges=graph.array_context(ctx.repeat,ctx.spacer);shared.insert(edges.first.begin(),edges.first.end());
    }
    for(const auto& cycle:cycles)source_models.push_back({"","",cycle});
    SourceCycleGuard guard(graph.k(),&source_models);
    auto candidates=std::move(result.alternatives);
    if(candidates.empty())candidates.push_back(result.selected);
    std::vector<Candidate> accepted,diagnostic;
    std::set<std::tuple<std::string,size_t,size_t>> seen;
    auto assess=[&](Candidate c) {
        size_t begin=c.seed_begin+result.placement.query_offset,end=begin+query.size();
        if(!seen.emplace(c.sequence,begin,end).second)return;
        auto q=assess_path_quality(graph,c,shared,guard,0);
        if(end<=c.sequence.size() && std::max(begin,c.sequence.size()-end)<options.min_extension_flank)
            q.reasons.push_back("requested_flank_not_reached");
        if(end>c.sequence.size() || c.sequence.substr(begin,query.size())!=query)
            q.reasons.push_back("target_sequence_differs");
        result.quality_checks.push_back({c.sequence,begin,end,q});
        c.score=q.score;
        if(q.accepted())accepted.push_back(c);
        if(options.report_alternatives)diagnostic.push_back(std::move(c));
    };
    for(auto& c:candidates) {
        // A bad flank cannot conceal a usable other side, nor can a good flank
        // make a bad concatenated region pass. Test the one-sided paths too.
        if(c.seed_begin && c.sequence.size()>c.seed_end) {
            auto left=c;trim_right(left,graph.k());left.right.stop="quality_single_flank";assess(std::move(left));
            auto right=c;trim_left(right);right.left.stop="quality_single_flank";assess(std::move(right));
        }
        assess(std::move(c));
    }
    auto better=[](const Candidate& a,const Candidate& b) {
        if(a.score!=b.score)return a.score>b.score;
        if(a.sequence.size()!=b.sequence.size())return a.sequence.size()>b.sequence.size();
        return a.sequence<b.sequence;
    };
    std::sort(accepted.begin(),accepted.end(),better);
    std::set<std::string> passing;for(const auto& c:accepted)passing.insert(c.sequence);
    std::stable_sort(diagnostic.begin(),diagnostic.end(),[&](const auto& a,const auto& b) {
        if(passing.count(a.sequence)!=passing.count(b.sequence))return passing.count(a.sequence)>passing.count(b.sequence);
        return better(a,b);
    });
    result.alternatives=std::move(diagnostic);
    result.has_runner=accepted.size()>1;result.runner_score=result.has_runner?accepted[1].score:0;
    result.has_extension=!accepted.empty();
    if(accepted.empty()) {
        result.selected={};result.status="not_extended";result.reason="no_path_passed_quality_filter";result.ambiguous=false;return;
    }
    result.selected=std::move(accepted.front());
    for(size_t i=1;i<accepted.size();++i)
        if(accepted[i].score>=result.selected.score-options.prediction_margin && distinct_extension_routes(result.selected,accepted[i],graph.k()))
            ++result.competing_quality_routes;
    if(result.competing_quality_routes) {
        result.quality_unresolved=true;result.has_extension=false;result.status="not_extended";
        result.reason="competing_quality_paths";result.ambiguous=true;return;
    }
    result.ambiguous=result.has_runner && result.runner_score>=result.selected.score-options.prediction_margin;
    result.status=result.variant_placement?(result.ambiguous?"extended_variant_ambiguous":"extended_variant"):
        (result.ambiguous?"extended_ambiguous":"extended");
    result.reason="quality_filter_pass";
}
} // namespace mcaat::proto
#endif
