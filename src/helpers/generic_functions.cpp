#include "generic_functions.h"
using namespace std;


/// @brief Helper Function gets all outgoing edges of a node
/// @param uint64_t node
string FetchNodeLabel(SDBG& succinct_de_bruijn_graph,auto node){
    std::string label;            
    uint8_t seq[succinct_de_bruijn_graph.k()];
    uint32_t t = succinct_de_bruijn_graph.GetLabel(node, seq);
    for (int i = succinct_de_bruijn_graph.k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

/// @brief Helper Function gets all outgoing edges of a node
/// @param uint64_t node
/// @return Outgoing edges of a node
set<uint64_t> GetOutgoingNodes(SDBG& succinct_de_bruijn_graph, uint64_t node){
    if (succinct_de_bruijn_graph.EdgeOutdegree(node) == 0 || node == -1)
        return {};
    int edge_outdegree = succinct_de_bruijn_graph.EdgeOutdegree(node);
    uint64_t outgoings[edge_outdegree];
    succinct_de_bruijn_graph.OutgoingEdges(node,outgoings);
    set<uint64_t> outgoings_set(outgoings, outgoings + edge_outdegree);
    return outgoings_set;
};
        
/// @brief Helper Function gets all incoming edges of a node
/// @param uint64_t node
/// @return Incoming edges of a node
set<uint64_t> GetIncomingNodes(SDBG& succinct_de_bruijn_graph,uint64_t node){
    if (succinct_de_bruijn_graph.EdgeIndegree(node) == 0 || node == -1)
        return {};
            
    int edge_indegree = succinct_de_bruijn_graph.EdgeIndegree(node);
    uint64_t incomings[edge_indegree];
    succinct_de_bruijn_graph.IncomingEdges(node,incomings);
    set<uint64_t> incomings_set;
    for (int i = 0; i < edge_indegree; i++){ 
        incomings_set.insert(incomings[i]);
    }
    
    return incomings_set;
};