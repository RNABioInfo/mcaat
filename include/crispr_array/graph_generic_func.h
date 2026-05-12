#ifndef GRAPH_GENERIC_FUNC_H
#define GRAPH_GENERIC_FUNC_H

#include <unordered_set>
#include <cstdint>
#include <sdbg/sdbg.h>
// Forward declaration for sdbg if needed
struct sdbg;

struct graph_generic_func
{
    //#### HELPER FUNCTIONS FOR DLS ###
    static void _GetOutgoings(uint64_t node, std::unordered_set<uint64_t>& outgoings_set, SDBG& sdbg);
    static void _GetIncomings(uint64_t node, std::unordered_set<uint64_t>& incomings_set, SDBG& sdbg);
    //#### HELPER FUNCTIONS FOR DLS ####
};

#endif // GRAPH_GENERIC_FUNC_H