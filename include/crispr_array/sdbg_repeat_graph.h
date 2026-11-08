#ifndef MCAAT_SDBG_REPEAT_GRAPH_H
#define MCAAT_SDBG_REPEAT_GRAPH_H
#include "crispr_array/graph_repeat_filter.h"
#include "sdbg/sdbg.h"
#include <sstream>

namespace mcaat::arrays {
// A k-mer vertex can have several outgoing SDBG edge slots. Normalize its
// identity, but retain the actual traversed edges separately for downstream code.
class SdbgRepeatGraph {
    const SDBG& g_;
    std::map<Node,std::string> labels_;
    std::map<Node,Node> reverse_;
public:
    explicit SdbgRepeatGraph(const SDBG& graph):g_(graph) {}
    size_t k() const {return g_.k();}
    std::string label(Node vertex) {
        auto found=labels_.find(vertex);
        if (found!=labels_.end()) return found->second;
        if (vertex>=g_.size()) throw std::runtime_error("graph vertex out of range");
        std::vector<uint8_t> bases(k());
        g_.GetLabel(vertex,bases.data());
        std::string result;
        for (auto base:bases) {
            if (base<1 || base>4) throw std::runtime_error("non-DNA graph label");
            result+="ACGT"[base-1];
        }
        labels_[vertex]=result;return result;
    }
    Node find(const std::string& sequence) const {
        if (sequence.size()!=k()) return missing_node;
        std::vector<uint8_t> bases;
        for (char c:sequence) {
            auto pos=std::string("ACGT").find(c);
            if (pos==std::string::npos) return missing_node;
            bases.push_back(uint8_t(pos+1));
        }
        return g_.IndexBinarySearch(bases.data());
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
