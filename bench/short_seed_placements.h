#ifndef MCAAT_SHORT_SEED_PLACEMENTS_H
#define MCAAT_SHORT_SEED_PLACEMENTS_H
#include "protospacer_detection/extension_engine.h"
namespace mcaat::proto::audit {
// Enumerate all k-mer contexts compatible with an immutable short input.
// The other side's context bases are not emitted as recovered flank sequence.
template<class Graph> Path short_seed_endpoints(Graph& graph,const std::string& seed,bool backward) {
    if(seed.empty() || seed.size()>=graph.k())throw std::runtime_error("expected short seed");
    size_t extra=graph.k()-seed.size();
    if(extra>8)throw std::runtime_error("short-seed context enumeration exceeds diagnostic limit");
    Path result;std::string context(extra,'A');
    for(size_t code=0;code<(size_t(1)<<(2*extra));++code) {
        auto value=code;for(auto& c:context){c="ACGT"[value&3];value>>=2;}
        auto word=backward?seed+context:context+seed;
        auto vertex=graph.find(word);if(vertex!=missing)result.push_back(vertex);
    }
    return result;
}
}
#endif
