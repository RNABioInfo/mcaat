#ifndef MCAAT_TARGET_PLACEMENT_H
#define MCAAT_TARGET_PLACEMENT_H
#include "protospacer_detection/extension_engine.h"
#include <functional>
namespace mcaat::proto {
struct TargetPlacement {
    std::string sequence; // actual graph sequence, including short-query padding
    size_t query_offset=0,query_begin=0,query_end=0,mismatches=0;
    double identity=0,coverage=0;
};
struct TargetPlacements {std::vector<TargetPlacement> values;size_t work=0;bool limited=false;};
inline TargetPlacement score_target_placement(const std::string& query,const std::string& observed,size_t offset,
                                               double min_identity=.90,double min_coverage=.85) {
    TargetPlacement result;result.sequence=observed;result.query_offset=offset;
    if(observed.size()<offset+query.size())return result;
    size_t minimum=size_t(std::ceil(query.size()*min_coverage-1e-9));
    for(size_t length=query.size();length>=minimum && length>0;--length) {
        bool found=false;
        for(size_t begin=0;begin+length<=query.size();++begin) {
            size_t errors=0;for(size_t i=begin;i<begin+length;++i)errors+=query[i]!=observed[offset+i];
            double identity=double(length-errors)/length;
            if(identity+1e-12<min_identity)continue;
            if(!found || errors<result.mismatches) {
                found=true;result.query_begin=begin;result.query_end=begin+length;result.mismatches=errors;
                result.identity=identity;result.coverage=double(length)/query.size();
            }
        }
        if(found)break;
    }
    return result;
}
// Graph-only, substitution-aware fallback. Locate exact/one-substitution k-mer
// anchors (or suffix anchors for queries shorter than k), then walk actual
// edges to cover the query. Final candidates must align >=85% of the query at
// >=90% identity. No reference bases are supplied or synthesized into a path.
// The bounded search is explicitly not a complete indel-aware mapper.
template<class Graph> TargetPlacements place_target(Graph& graph,const std::string& query,size_t left_padding,
                                                    size_t max_work=50000,size_t max_candidates=128) {
    TargetPlacements result;if(query.size()<15)return result;
    const size_t k=graph.k(),core=std::max(k,query.size());
    if(left_padding>core-query.size())throw std::runtime_error("invalid target padding");
    const size_t allowed=query.size()-size_t(std::ceil(.85*query.size()-1e-9))+size_t(std::floor(.10*query.size()+1e-9));
    std::set<std::string> seen;std::set<std::pair<Node,size_t>> anchors;
    auto work=[&](){if(result.work>=max_work){result.limited=true;return false;}++result.work;return true;};
    auto accept=[&](const Path& path) {
        auto sequence=graph.spell(path);if(sequence.size()!=core)throw std::runtime_error("target placement length mismatch");
        if(!seen.insert(sequence).second)return;
        auto scored=score_target_placement(query,sequence,left_padding);
        if(!scored.query_end)return;
        if(result.values.size()>=max_candidates){result.limited=true;return;}
        result.values.push_back(std::move(scored));
    };
    auto errors=[&](Node node,int64_t begin) {
        size_t n=0;auto word=graph.label(node);
        for(size_t j=0;j<word.size();++j)if(begin+int64_t(j)>=0 && begin+int64_t(j)<int64_t(query.size()))
            n+=word[j]!=query[size_t(begin+j)];
        return n;
    };
    auto expand=[&](Node anchor,size_t anchor_pos) {
        if(!anchors.emplace(anchor,anchor_pos).second)return;
        if(query.size()<k) {
            // The suffix anchor ends at query position anchor_pos. Shift to a
            // complete observed k-mer around the query; every padding base is
            // read from a real edge and is reported outside the matched input.
            size_t steps=k-left_padding-anchor_pos;
            std::function<void(Node,size_t)> walk=[&](Node node,size_t remaining) {
                if(!work() || errors(node,-int64_t(left_padding)-int64_t(remaining))>allowed)return;
                if(!remaining){accept(Path{node});return;}
                for(auto arc:graph.neighbors(node,false))walk(arc.to,remaining-1);
            };
            walk(anchor,steps);return;
        }
        Path path{anchor};
        std::function<void(Path,size_t,size_t)> right=[&](Path nodes,size_t position,size_t mismatches) {
            if(!work() || mismatches>allowed)return;
            if(position==query.size()){accept(nodes);return;}
            for(auto arc:graph.neighbors(nodes.back(),false)) {
                auto next=nodes;next.push_back(arc.to);
                right(std::move(next),position+1,mismatches+(graph.extension_base(arc,false)!=query[position]));
            }
        };
        std::function<void(Path,size_t,size_t)> left=[&](Path nodes,size_t position,size_t mismatches) {
            if(!work() || mismatches>allowed)return;
            if(!position){right(std::move(nodes),anchor_pos+k,mismatches);return;}
            for(auto arc:graph.neighbors(nodes.front(),true)) {
                auto next=nodes;next.insert(next.begin(),arc.to);
                left(std::move(next),position-1,mismatches+(graph.extension_base(arc,true)!=query[position-1]));
            }
        };
        left(path,anchor_pos,errors(anchor,int64_t(anchor_pos)));
    };
    const size_t anchor_length=query.size()<k?size_t(std::ceil(.85*query.size()-1e-9)):k;
    // Examine exact anchors before substitutions, across all query positions.
    for(size_t edits=0;edits<=1 && !result.limited;++edits)
      for(size_t position=0;position+anchor_length<=query.size() && !result.limited;++position) {
        auto anchor=query.substr(position,anchor_length);
        auto lookup=[&](const std::string& word) {
            if(!work())return;
            if(query.size()<k) {
                bool capped=false;auto nodes=graph.find_suffix(word,max_candidates,capped);result.limited|=capped;
                for(auto node:nodes)expand(node,position+anchor_length);
            }else {auto node=graph.find(word);if(node!=missing)expand(node,position);}
        };
        if(!edits)lookup(anchor);
        else for(size_t i=0;i<anchor.size()&&!result.limited;++i)for(char base:std::string("ACGT"))if(base!=anchor[i]) {
            auto word=anchor;word[i]=base;lookup(word);if(result.limited)break;
        }
      }
    std::sort(result.values.begin(),result.values.end(),[](const auto& a,const auto& b){
        if(a.coverage!=b.coverage)return a.coverage>b.coverage;
        if(a.identity!=b.identity)return a.identity>b.identity;
        return a.sequence<b.sequence;
    });
    return result;
}
}
#endif
