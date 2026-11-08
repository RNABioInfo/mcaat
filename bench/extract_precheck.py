"""Wrap original MCAAT methods for benchmarking against the actual SDBG backend."""
import pathlib
import sys

source = pathlib.Path(sys.argv[1]).read_text()
names = [
    "bool CycleFinder::BidirectionalBFS(uint64_t start, int limit)",
    "vector<uint64_t> CycleFinder::CollectTips()",
    "void CycleFinder::RecursiveReduction(uint64_t tip)",
    "void CycleFinder::InvalidateMultiplicityOneNodes()",
    "void CycleFinder::_GetIncomings(uint64_t node, SmallNodeSet& incomings_set)",
]
methods = []
for signature in names:
    begin = source.index(signature)
    brace = source.index("{", begin)
    depth = 1
    end = brace + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    methods.append(source[begin:end].replace("CycleFinder::", "OriginalPrecheck::"))

wrapper = r'''
#include <stack>
#include <unordered_set>
#include <limits>
#include <omp.h>
using std::vector;
struct SmallNodeSet {
    uint64_t nodes[4]; uint8_t count = 0;
    void insert(uint64_t v) { nodes[count++] = v; }
    uint64_t* begin() { return nodes; }
    uint64_t* end() { return nodes + count; }
};
class OriginalPrecheck {
    struct Settings { SDBG* sdbg; size_t threads = 1; } settings;
    static constexpr int MAX_EDGE_COUNT = 4;
    vector<vector<uint64_t>> per_thread_visited;
    vector<vector<uint32_t>> per_thread_dirty;
    void _GetIncomings(uint64_t node, SmallNodeSet& incomings_set);
public:
    explicit OriginalPrecheck(SDBG& g) : settings{&g},
        per_thread_visited(1,vector<uint64_t>((g.size()+63)/64)),per_thread_dirty(1) {}
    bool BidirectionalBFS(uint64_t start, int limit);
    vector<uint64_t> CollectTips();
    void RecursiveReduction(uint64_t tip);
    void InvalidateMultiplicityOneNodes();
};
'''
pathlib.Path(sys.argv[2]).write_text(wrapper + "\n\n".join(methods) + "\n")
