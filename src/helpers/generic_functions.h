#ifndef GENERIC_FUNCTIONS_H
#define GENERIC_FUNCTIONS_H

#include <vector>
#include <map>
#include <set>
#include <array>
#include <../src/sdbg/sdbg.h>
#include <memory>
using namespace std;

string FetchNodeLabel(SDBG& succinct_de_bruijn_graph,auto node);

set<uint64_t> GetOutgoingNodes(SDBG& succinct_de_bruijn_graph, uint64_t node);

set<uint64_t> GetIncomingNodes(SDBG& succinct_de_bruijn_graph,uint64_t node);
#endif