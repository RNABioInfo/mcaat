#ifndef MCAAT_COMPLETE_TARGET_PLACEMENTS_H
#define MCAAT_COMPLETE_TARGET_PLACEMENTS_H
#include "protospacer_detection/target_placement.h"

namespace mcaat::proto::audit {
struct CompletePlacement {Node endpoint;std::string observed,context,flank_prefix;size_t begin,end,errors;};
struct CompletePlacements {std::map<std::pair<Node,std::string>,CompletePlacement> endpoints;size_t lookups=0,states=0;};
// Exhaustive for the production mapper's UNGAPPED >=90% identity / >=85%
// query-coverage policy. No candidate, work, or multiplicity cutoff. Search all
// graph k-mers satisfying each permitted query interval, then validate all
// connected continuations. Short inputs can occur anywhere inside a k-mer;
// flank_prefix records the already observed bases after such an input.
template<class Graph> CompletePlacements complete_target_placements(Graph& graph,const std::string& query) {
    CompletePlacements result;const size_t n=query.size(),k=graph.k();
    if(!n)throw std::runtime_error("empty query");
    std::set<size_t> lengths{(85*n+99)/100};
    // Any longer qualifying interval contains an interval of one of these
    // lengths with the same substitution budget. No valid endpoint is lost.
    for(size_t length=*lengths.begin();length<=n;++length)if(length%10==0)lengths.insert(length);
    // A single relaxed index lookup covers all allowed query intervals. The
    // common interval is present in every >=85% alignment; all other positions
    // are wildcards. Full interval checks below remove its extra candidates.
    std::vector<std::vector<std::pair<Node,size_t>>> candidates;
    for(size_t offset=0;offset<=(n<k?k-n:0);++offset) {
        std::string pattern(k,'.');size_t minimum=*lengths.begin(),clip=n-minimum;
        for(size_t j=offset;j<k&&j-offset<n;++j)if(j-offset>=clip&&j-offset<minimum)pattern[j]=query[j-offset];
        ++result.lookups;candidates.push_back(graph.find_masked(pattern,n/10));
    }
    for(auto length:lengths)for(size_t begin=0;begin+length<=n;++begin) {
        size_t end=begin+length,allowed=length/10;
        for(size_t offset=0;offset<=(n<k?k-n:0);++offset) {
            const auto& roots=candidates[offset];
            struct Cell {size_t errors;Path path;};
            std::map<Node,Cell> current;
            for(auto [node,ignored]:roots) {
                auto label=graph.label(node);size_t errors=0;
                for(size_t j=offset;j<k&&j-offset<n;++j)if(j-offset>=begin&&j-offset<end)errors+=label[j]!=query[j-offset];
                if(errors<=allowed)current.emplace(node,Cell{errors,Path{node}});
            }
            for(size_t pos=k;pos<n&&!current.empty();++pos) {
                std::map<Node,Cell> next;
                for(const auto& [v,cell]:current)for(auto arc:graph.neighbors(v,false)) {
                    ++result.states;size_t cost=cell.errors+(pos>=begin&&pos<end?graph.extension_base(arc,false)!=query[pos]:0);
                    if(cost>allowed)continue;
                    auto found=next.find(arc.to);if(found!=next.end()&&found->second.errors<=cost)continue;
                    auto path=cell.path;path.push_back(arc.to);next[arc.to]={cost,std::move(path)};
                }
                current=std::move(next);
            }
            for(const auto& [v,cell]:current) {
                auto context=graph.spell(cell.path),observed=context.substr(offset,n),flank=context.substr(offset+n);
                auto score=score_target_placement(query,observed,0);
                if(!score.query_end)throw std::runtime_error("complete mapper emitted invalid query placement");
                auto key=std::make_pair(v,flank);auto found=result.endpoints.find(key);
                if(found==result.endpoints.end() || score.coverage>double(found->second.end-found->second.begin)/n ||
                        (score.query_end-score.query_begin==found->second.end-found->second.begin&&score.mismatches<found->second.errors))
                    result.endpoints[key]={v,observed,context,flank,score.query_begin,score.query_end,score.mismatches};
            }
        }
    }
    return result;
}
}
#endif
