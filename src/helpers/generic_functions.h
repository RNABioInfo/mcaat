#ifndef GENERIC_FUNCTIONS_H
#define GENERIC_FUNCTIONS_H

#include <vector>
#include <map>
#include <set>
#include <array>
#include <../src/sdbg/sdbg.h>

using namespace std;

string FetchNodeLabel(SDBG& succinct_de_bruijn_graph,auto node);

set<uint64_t> GetOutgoingNodes(SDBG& succinct_de_bruijn_graph, uint64_t node);

set<uint64_t> GetIncomingNodes(SDBG& succinct_de_bruijn_graph,uint64_t node);

void GenericWrite(vector<uint64_t> path, uint64_t startnode, string type, string genome_id, 
    SDBG& succinct_de_bruijn_graph, char file_mode,string folder_path="proof_of_concept/data/");

void WritePath(vector<uint64_t> path, uint64_t startnode, string modes, string type,
 string genome_id, SDBG& succinct_de_bruijn_graph,string folder_path="proof_of_concept/data/");

void WriteBlankLine(string type,string genome_id,string folder_path="proof_of_concept/data/");

#endif