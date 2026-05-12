#include <unordered_set>
#include "crispr_array/graph_generic_func.h"


// generic helper functions for graph operations

void graph_generic_func::_GetOutgoings(uint64_t node, std::unordered_set<uint64_t>& outgoings_set, SDBG& sdbg) {
    int edge_outdegree = sdbg.EdgeOutdegree(node);
    if (edge_outdegree == 0 || !sdbg.IsValidEdge(node)) {
        return;
    }
    std::vector<uint64_t> outgoings(edge_outdegree);
    int flag = sdbg.OutgoingEdges(node, outgoings.data());
    if (flag != -1)
        for (const auto& outgoing : outgoings)
            if(sdbg.EdgeMultiplicity(outgoing) >= sdbg.EdgeMultiplicity(node)/2 && sdbg.EdgeMultiplicity(outgoing) <= sdbg.EdgeMultiplicity(node)*1.2 && sdbg.IsValidEdge(outgoing))
                // filter out edges with higher variability
                outgoings_set.insert(outgoing);
}

void graph_generic_func::_GetIncomings(uint64_t node, std::unordered_set<uint64_t>& incomings_set, SDBG& sdbg) {
    int edge_indegree = sdbg.EdgeIndegree(node);
    if (edge_indegree == 0 || !sdbg.IsValidEdge(node)) {
        return;
    }
    std::vector<uint64_t> incomings(edge_indegree);
    int flag = sdbg.IncomingEdges(node, incomings.data());
    if (flag != -1)
        for (const auto& incoming : incomings)
            if(sdbg.EdgeMultiplicity(incoming) >= sdbg.EdgeMultiplicity(node)/2 && 
            sdbg.EdgeMultiplicity(incoming) <= sdbg.EdgeMultiplicity(node)*1.2 &&
             sdbg.IsValidEdge(incoming))
                // filter out edges with higher variability
                incomings_set.insert(incoming);
            
}