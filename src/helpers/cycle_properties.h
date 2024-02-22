//generate a structure to hold the properties of a graph
//this is a helper structure for the graph class

#ifndef GRAPH_PROPERTIES_H
#define GRAPH_PROPERTIES_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <../src/sdbg/sdbg.h>

using namespace std;

struct CycleProperties { 
        uint8_t maximal_length;
        uint8_t minimal_length;
        SDBG& succinct_de_bruijn_graph;
        string genome_name;
        map<uint64_t, bool> visited;
        uint8_t cluster_bounds;
};

#endif

