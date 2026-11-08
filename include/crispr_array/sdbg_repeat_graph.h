#ifndef MCAAT_SDBG_REPEAT_GRAPH_H
#define MCAAT_SDBG_REPEAT_GRAPH_H
#include "crispr_array/graph_repeat_filter.h"
#include "sdbg/sdbg.h"
#include <sstream>
#include <unordered_map>

namespace mcaat::arrays {
// A k-mer vertex can have several outgoing SDBG edge slots. Normalize its
// identity, but retain the actual traversed edges separately for downstream code.
class SdbgRepeatGraph {
    const SDBG& g_;
    std::unordered_map<Node,std::string> labels_;
    std::unordered_map<Node,Node> reverse_;
    std::unordered_map<std::string,Node> lookups_;
    static constexpr size_t lookup_cache_limit=65536;
public:
    size_t index_queries=0,lookup_hits=0,walk_steps=0;
    explicit SdbgRepeatGraph(const SDBG& graph):g_(graph) {}
    size_t k() const {return g_.k();}
    const std::string& label(Node vertex) {
        auto found=labels_.find(vertex);
        if (found!=labels_.end()) return found->second;
        if (vertex>=g_.size()) throw std::runtime_error("graph vertex out of range");
        std::vector<uint8_t> bases(k());
        g_.GetLabel(vertex,bases.data());
        std::string result;
        result.reserve(k());
        for (auto base:bases) {
            if (base<1 || base>4) throw std::runtime_error("non-DNA graph label");
            result+="ACGT"[base-1];
        }
        return labels_.emplace(vertex,std::move(result)).first->second;
    }
    Node find(const std::string& sequence) {
        if (sequence.size()!=k()) return missing_node;
        auto cached=lookups_.find(sequence);
        if (cached!=lookups_.end()) {++lookup_hits;return cached->second;}
        std::vector<uint8_t> bases;
        bases.reserve(k());
        for (char c:sequence) {
            switch(c) {
                case 'A':bases.push_back(1);break;case 'C':bases.push_back(2);break;
                case 'G':bases.push_back(3);break;case 'T':bases.push_back(4);break;
                default:return missing_node;
            }
        }
        ++index_queries;
        Node vertex=g_.IndexBinarySearch(bases.data());
        // The start/RC lookup cache is bounded independently of input size.
        if (lookups_.size()>=lookup_cache_limit) lookups_.clear();
        lookups_.emplace(sequence,vertex);
        return vertex;
    }
    // A verified suffix/prefix overlap identifies the next k-mer by its final
    // base. Follow that edge instead of searching the entire graph again.
    std::pair<Node,Node> step(Node from,char base) {
        ++walk_steps;
        uint8_t encoded;
        switch(base) {
            case 'A':encoded=1;break;case 'C':encoded=2;break;
            case 'G':encoded=3;break;case 'T':encoded=4;break;
            default:return {missing_node,missing_node};
        }
        if (from>=g_.size()) return {missing_node,missing_node};
        Node slot=from;
        do {
            auto w=g_.GetW(slot);
            if (w>4) w-=4;
            if (w==encoded && g_.IsValidEdge(slot)) return {slot,g_.Forward(slot)};
            --slot;
        } while (slot!=missing_node && !g_.IsLastOrTip(slot));
        return {missing_node,missing_node};
    }
    Node edge(Node from,Node to) const {
        if (from>=g_.size() || to>=g_.size()) return missing_node;
        Node slot=from;
        do {
            if (g_.IsValidEdge(slot) && g_.Forward(slot)==to) return slot;
            --slot;
        } while (slot!=missing_node && !g_.IsLastOrTip(slot));
        return missing_node;
    }
    bool adjacent(Node from,Node to) const {return edge(from,to)!=missing_node;}
    Node reverse_node(Node vertex) {
        auto found=reverse_.find(vertex);
        if (found!=reverse_.end()) return found->second;
        const Node rc=find(reverse_complement(label(vertex)));
        reverse_[vertex]=rc;return rc;
    }
    Path vertices(const Path& edges) const {
        Path result;
        result.reserve(edges.size());
        for (size_t i=0;i<edges.size();++i) {
            Node e=edges[i],next=edges[(i+1)%edges.size()];
            if (e>=g_.size() || next>=g_.size() || !g_.IsValidEdge(e) ||
                g_.Forward(e)!=uint64_t(g_.GetLastIndex(next))) return {};
            result.push_back(g_.GetLastIndex(e));
        }
        return result;
    }
    Path edges(const Path& vertices) const {
        Path result;
        result.reserve(vertices.size());
        for (size_t i=0;i<vertices.size();++i) {
            Node e=edge(vertices[i],vertices[(i+1)%vertices.size()]);
            if (e==missing_node) throw std::runtime_error("cycle is not present in graph");
            result.push_back(e);
        }
        return result;
    }
};

inline std::string path_text(const Path& path) {
    if (path.empty()) return "-";
    std::ostringstream out;
    for (size_t i=0;i<path.size();++i) {if(i) out<<',';out<<path[i];}
    return out.str();
}
inline Path parse_path(const std::string& text) {
    if (text=="-") return {};
    Path result;std::istringstream input(text);std::string token;
    while (std::getline(input,token,',')) {
        if (token.empty() || token.find_first_not_of("0123456789")!=std::string::npos)
            throw std::runtime_error("invalid graph path ID");
        result.push_back(std::stoull(token));
    }
    if (result.empty() || text.back()==',') throw std::runtime_error("empty graph path");
    return result;
}
} // namespace mcaat::arrays
#endif
