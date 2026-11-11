#!/usr/bin/env python3
"""Copy current production headers and inject logging only, outside the repository."""
import argparse
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('output',type=Path);a=p.parse_args()
repo=Path(__file__).resolve().parent.parent
include=a.output/'instrumented/protospacer_detection';include.mkdir(parents=True,exist_ok=True)
s=(repo/'include/protospacer_detection/extension_engine.h').read_text()
def replace(old,new):
 global s
 assert s.count(old)==1,(old[:100],s.count(old))
 s=s.replace(old,new)
replace('    auto finish=[&](Side state,std::string reason) {','    auto finish=[&](Side state,std::string reason) {\n        audit_side("stopping_prefix",backward,state,reason);')
replace('        state.stop=std::move(reason);','        audit_side("terminal",backward,state,reason);\n        state.stop=std::move(reason);')
replace('if(!stop.empty()){repeat_reason=stop;continue;}','if(!stop.empty()){audit_arc("reject",backward,current,a,0,stop);repeat_reason=stop;continue;}')
replace('if(array_edges.count(a.edge)){masked=true;continue;}','if(array_edges.count(a.edge)){audit_arc("reject",backward,current,a,0,"masked_array_edge");masked=true;continue;}')
replace('                if(!adapter.empty()) {\n                    // The completing base', '                if(!adapter.empty()) {\n                    audit_arc("reject",backward,next,a,0,"sequencing_adapter");\n                    // The completing base')
replace('                        next.context="array";rejected.emplace_back(std::move(next),"source_array_context");continue;', '                        audit_arc("reject",backward,next,a,0,"source_array_context");\n                        next.context="array";rejected.emplace_back(std::move(next),"source_array_context");continue;')
replace('if(penalty>options.max_coverage_penalty) {next.coverage=std::move(before);rejected.emplace_back(std::move(next),"coverage_break");continue;}', 'if(penalty>options.max_coverage_penalty) {next.coverage=std::move(before);audit_arc("reject",backward,next,a,penalty,"coverage_break");rejected.emplace_back(std::move(next),"coverage_break");continue;}')
replace('            next_frontier.resize(options.beam_width);result.beam_pruned=true;', '            for(size_t z=options.beam_width;z<next_frontier.size();++z)audit_side("beam_pruned",backward,next_frontier[z],std::to_string(z));\n            next_frontier.resize(options.beam_width);result.beam_pruned=true;')
replace('if(result.paths.size()>options.beam_width){result.paths.resize(options.beam_width);result.beam_pruned=true;}', 'if(result.paths.size()>options.beam_width){for(size_t z=options.beam_width;z<result.paths.size();++z)audit_side("terminal_pruned",backward,result.paths[z],result.paths[z].stop);result.paths.resize(options.beam_width);result.beam_pruned=true;}')
replace('    if(unique.size()>options.beam_width){unique.resize(options.beam_width);result.beam_pruned=true;}', '    if(unique.size()>options.beam_width){for(size_t z=options.beam_width;z<unique.size();++z)audit_side("side_shortlist_pruned",backward,unique[z],unique[z].stop);unique.resize(options.beam_width);result.beam_pruned=true;}')
replace('    result.paths=std::move(unique);return result;', '    result.paths=std::move(unique);for(auto& side:result.paths)audit_side("side_retained",backward,side,side.stop);return result;')
replace('if(!compatible_walk(graph,c.path,seed.path,ln,join_guard,options))continue;', 'if(!compatible_walk(graph,c.path,seed.path,ln,join_guard,options)){audit_candidate("join_rejected",c,0);continue;}')
replace('        results.push_back(std::move(c));\n    }\n    if(results.size()>options.beam_width)', '        audit_candidate("joined",c,0);results.push_back(std::move(c));\n    }\n    if(results.size()>options.beam_width)')
replace('        results.resize(options.beam_width);for(auto& c:results)c.search_limited=true;', '        for(size_t z=options.beam_width;z<results.size();++z)audit_candidate("join_shortlist_pruned",results[z],z);\n        results.resize(options.beam_width);for(auto& c:results)c.search_limited=true;')
(include/'extension_engine.h').write_text(s)
s=(repo/'include/protospacer_detection/retained_target.h').read_text()
replace('        refresh_target(graph,c,left_models,right_models,shared,options);', '        audit_candidate("before_refresh",c,0);\n        refresh_target(graph,c,left_models,right_models,shared,options);\n        audit_candidate("after_refresh",c,0);')
replace('            result.stops.insert(refined.left.stop);result.stops.insert(refined.right.stop);','            audit_candidate("after_refine",refined,0);\n            result.stops.insert(refined.left.stop);result.stops.insert(refined.right.stop);')
replace('    const auto& winner=candidates.front();','    for(size_t rank=0;rank<candidates.size();++rank)audit_candidate("eligible",candidates[rank],rank);\n    const auto& winner=candidates.front();')
(include/'retained_target.h').write_text(s)
print('Logging-only headers:',include)
