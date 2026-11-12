#ifndef MCAAT_REFERENCE_SEGMENT_BOUNDS_H
#define MCAAT_REFERENCE_SEGMENT_BOUNDS_H
#include "protospacer_detection/extension_engine.h"

namespace mcaat::proto::audit {
// A lower bound, not a mapper: allow any graph vertex and even ignore the
// connecting edge between the two k-mers of a k+1-base word. This relaxation
// can only make the bound smaller than the true alignment cost.
template<class Graph> size_t segment_lower_bound(Graph& graph,const std::string& word) {
    const size_t k=graph.k();
    if(word.size()!=k)throw std::runtime_error("segment bound requires k bases");
    if(graph.find(word)!=missing)return 0;
    // One insertion costs 9. Both end k-mers must exist; their connecting edge
    // need not exist in this deliberately optimistic membership test.
    for(size_t i=0;i<=k;++i)for(char base:std::string("ACGT")) {
        auto candidate=word.substr(0,i)+base+word.substr(i);
        if(graph.find(candidate.substr(0,k))!=missing && graph.find(candidate.substr(1))!=missing)return 9;
    }
    // One substitution or deletion costs 10. A k-1-base emitted segment must
    // be the suffix of some reached k-mer, including during a backwards walk.
    for(size_t i=0;i<k;++i) {
        auto candidate=word;
        for(char base:std::string("ACGT"))if(base!=word[i]) {
            candidate[i]=base;if(graph.find(candidate)!=missing)return 10;
        }
        candidate=word.substr(0,i)+word.substr(i+1);
        for(char base:std::string("ACGT"))if(graph.find(base+candidate)!=missing || graph.find(candidate+base)!=missing)return 10;
    }
    // Every alignment needs at least two edits, each costing at least 9.
    return 18;
}

// Stronger fallback for difficult cases. A k+2-base reference block aligned
// with <=2 edits emits at least k bases. Its first and last k-mers must exist.
// Ignore their internal connection, again making this an optimistic bound.
template<class Graph> size_t two_edit_segment_lower_bound(Graph& graph,const std::string& word) {
    const size_t k=graph.k();if(word.size()!=k+2)throw std::runtime_error("two-edit block requires k+2 bases");
    std::unordered_map<std::string,bool> membership;
    auto exists=[&](const std::string& s){auto it=membership.find(s);if(it==membership.end())it=membership.emplace(s,graph.find(s)!=missing).first;return it->second;};
    auto possible=[&](const std::string& s){return s.size()>=k && exists(s.substr(0,k)) && exists(s.substr(s.size()-k));};
    if(possible(word))return 0;
    auto edits=[&](const std::string& input,const auto& visit) {
        for(size_t i=0;i<=input.size();++i)for(char c:std::string("ACGT"))if(visit(input.substr(0,i)+c+input.substr(i),size_t(9)))return true;
        for(size_t i=0;i<input.size();++i) {
            if(visit(input.substr(0,i)+input.substr(i+1),size_t(10)))return true;
            auto s=input;for(char c:std::string("ACGT"))if(c!=input[i]){s[i]=c;if(visit(s,size_t(10)))return true;}
        }
        return false;
    };
    size_t cost=0;
    if(edits(word,[&](const std::string& s,size_t c){if(!possible(s))return false;cost=c;return true;}))return cost;
    if(edits(word,[&](const std::string& s,size_t){return edits(s,[&](const std::string& t,size_t){return possible(t);});}))return 18;
    return 27;
}

// Sum independent, disjoint reference blocks. Insertions at a block boundary
// are assigned to one adjacent block, so no edit is charged twice. Ignore a
// partial block at the current position and at the checkpoint.
struct SegmentBounds {
    size_t block=0;
    std::vector<size_t> costs;
    std::vector<size_t> prefix;
    void prepare(){prefix.assign(1,0);for(auto cost:costs)prefix.push_back(prefix.back()+cost);}
    size_t operator()(size_t pos,size_t target)const {
        size_t first=(pos+block-1)/block,last=target/block,total=0;
        if(first>=last)return 0;
        if(!prefix.empty())return prefix.at(last)-prefix.at(first);
        for(size_t i=first;i<last;++i)total+=costs.at(i);
        return total;
    }
};
} // namespace mcaat::proto::audit
#endif
