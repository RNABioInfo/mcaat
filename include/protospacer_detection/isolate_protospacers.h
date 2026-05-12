#ifndef ISOLATE_PROTOSPACERS_H
#define ISOLATE_PROTOSPACERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cstdint>
#include <functional>
#include "sdbg/sdbg.h"

//define pair<std::map<uint64_t,std::set<uint64_t>>,std::map<uint64_t,std::set<uint64_t>>>
typedef std::pair<std::map<uint64_t,std::set<uint64_t>>,std::map<uint64_t,std::set<uint64_t>>> IN_OUT_PAIR_MAP_SET;

class IsolateProtospacers {
private:
    const SDBG& graph;
    std::map<uint64_t, std::set<uint64_t>> cycleNodes;
    std::map<uint64_t, uint64_t> cycleToGroup;
    IN_OUT_PAIR_MAP_SET protospacerNodes;

    std::pair<std::map<uint64_t, std::set<uint64_t>>, std::map<uint64_t, std::set<uint64_t>>> FilterSharedGroups(const std::map<uint64_t, std::set<uint64_t>>& inGroup, const std::map<uint64_t, std::set<uint64_t>>& outGroup);

    void DepthLimitedSearch(uint64_t currentNode, int depth, std::vector<uint64_t>& path, std::set<uint64_t>& visited, const std::set<uint64_t>& outNodes, const std::set<uint64_t>& cycleNodeSet, int maxDepth, int minDepth, std::vector<std::vector<uint64_t>>& successfulPaths);

public:
    IsolateProtospacers(const SDBG& g, const std::string& cyclesFilePath);
    IsolateProtospacers(const SDBG& g, const std::map<uint64_t, std::vector<std::vector<uint64_t>>>& repeatToSpacerNodes);

    std::map<uint64_t, std::set<uint64_t>> ReadCycles(const std::string& filePath);

    IN_OUT_PAIR_MAP_SET IsolateProtospacersMethod();

    void PrintProtospacerNodesToConsole(const std::map<uint64_t, std::set<uint64_t>>& possibleProtospacerNodes, const std::string& type_edge);

    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> DepthLimitedPathsFromInToOut(
        const std::map<uint64_t, std::set<uint64_t>>& inGroup,
        const std::map<uint64_t, std::set<uint64_t>>& outGroup,
        int maxDepth,
        int minDepth
    );

    void WritePathsToFile(const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& paths, const std::string& filename);

    // Getters
    const IN_OUT_PAIR_MAP_SET& getProtospacerNodes() const { return protospacerNodes; }
    const std::map<uint64_t, std::set<uint64_t>>& getCycleNodes() const { return cycleNodes; }

    
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> ReadPathsFromFile(const std::string& filename);
};

#endif