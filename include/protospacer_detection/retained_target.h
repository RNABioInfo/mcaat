#ifndef MCAAT_PROTO_RETAINED_TARGET_H
#define MCAAT_PROTO_RETAINED_TARGET_H
#include "protospacer_detection/extension_engine.h"

namespace mcaat::proto {
// All strings are oriented like the supplied, immutable protospacer.
struct TargetContext {
    std::string spacer,repeat;
    size_t begin=0,end=0;
};
struct TargetExtension {
    std::string status="not_extended",reason;
    Candidate selected;
    bool has_extension=false,ambiguous=false,search_limited=false;
    size_t compared=0,expanded=0;
    double runner_score=0;
    bool has_runner=false;
    SearchDiagnostics diagnostics;
    std::set<std::string> stops;
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
        if(state!=ContextState::external || seq.size()<minimum) {
            if(backward)trim_left(c);else trim_right(c,graph.k());
            auto& side=backward?c.left:c.right;
            side.context=state==ContextState::array?"array":"incomplete";
            side.stop=state==ContextState::array?"source_array_context":
                state==ContextState::pending?"incomplete_context":"insufficient_external_sequence";
        }
    };
    check_side(left,left_models,true);check_side(right,right_models,false);
    c.external_left=c.external_right=0;
    for(size_t i=1;i<c.path.size();++i) {
        if(shared.count(graph.edge(c.path[i-1],c.path[i])))continue;
        if(i<=c.seed_begin)++c.external_left;
        if(i>c.seed_end-graph.k())++c.external_right;
    }
    if(!c.external_left && c.left_length){trim_left(c);c.left.stop="array_context_only";}
    if(!c.external_right && c.right_length){trim_right(c,graph.k());c.right.stop="array_context_only";}
    Path l(c.path.begin(),c.path.begin()+c.seed_begin+1);std::reverse(l.begin(),l.end());
    Path r(c.path.begin()+c.seed_end-graph.k(),c.path.end());
    c.left.coverage=path_coverage(graph,l,l.size()-1,true,c.left.context_clear_at,&shared,options);
    c.right.coverage=path_coverage(graph,r,r.size()-1,false,c.right.context_clear_at,&shared,options);
    c.score=stability(c.left.coverage,c.right.coverage).score();
}
template<class Graph>
TargetExtension extend_retained_target(Graph& graph,const std::string& sequence,
                                      const std::vector<TargetContext>& contexts,const Options& options,
                                      const std::set<std::string>* source_cycles=nullptr) {
    TargetExtension result;
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
        if(ctx.end<ctx.begin || ctx.end>ctx.spacer.size() || ctx.spacer.substr(ctx.begin,ctx.end-ctx.begin)!=sequence)
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
} // namespace mcaat::proto
#endif
