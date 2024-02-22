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
    if (succinct_de_bruijn_graph.EdgeOutdegree(node) == 0)
        return {};
            
    uint64_t *outgoings = new uint64_t[succinct_de_bruijn_graph.EdgeOutdegree(node)];
    succinct_de_bruijn_graph.OutgoingEdges(node,outgoings);
    set<uint64_t> outgoings_set(outgoings, outgoings + succinct_de_bruijn_graph.EdgeOutdegree(node));
    return outgoings_set;
};
        
/// @brief Helper Function gets all incoming edges of a node
/// @param uint64_t node
/// @return Incoming edges of a node
set<uint64_t> GetIncomingNodes(SDBG& succinct_de_bruijn_graph,uint64_t node){
    if (succinct_de_bruijn_graph.EdgeIndegree(node) == 0)
        return {};
            
    uint64_t *incomings = new uint64_t[succinct_de_bruijn_graph.EdgeIndegree(node)];
    succinct_de_bruijn_graph.IncomingEdges(node,incomings);
    set<uint64_t> incomings_set(incomings, incomings + succinct_de_bruijn_graph.EdgeIndegree(node));
    return incomings_set;
};

void Merge(int arr[], int ids[],int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];

    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];

    i = 0;
    j = 0;
    k = l;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void MergeSort(int arr[], int ids[], int n)
{
    int curr_size;
    int left_start;

    for (curr_size = 1; curr_size <= n - 1; curr_size = 2 * curr_size) {
        for (left_start = 0; left_start < n - 1; left_start += 2 * curr_size) {
            int mid = min(left_start + curr_size - 1, n - 1);
            int right_end = min(left_start + 2 * curr_size - 1, n - 1);

            Merge(arr, ids,left_start, mid, right_end);
        }
    }
}
