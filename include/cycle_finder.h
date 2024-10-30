#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H
#include <sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <array>
#include <stack>
#include <omp.h>
#include "path_writer.h"
#include <unordered_set>
//#include "progressbar.hpp"
#include <list>
#include <malloc.h>
using namespace std;

class CycleFinder {
    private:
        uint8_t maximal_length;
        uint8_t minimal_length;
        
        SDBG& sdbg;
        string genome_name;
        uint16_t cluster_bounds;
        vector<bool> visited;
        bool _IncomingNotEqualToCurrentNode(uint64_t node, size_t edge_indegree);
        bool _BackgroundCheck(auto original_node, size_t repeat_multiplicity, auto current_node,bool& stop);
        void _GetOutgoings(auto node, unordered_set<uint64_t>& outgoings_set, size_t repeat_multiplicity,bool& stop);
        void _GetIncomings(auto node, unordered_set<uint64_t>& incomings_set, size_t repeat_multiplicity,bool& stop);
        void _NeighborOut(auto node, unordered_set<uint64_t>& outgoings_set);
    public:
        CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name);
        
        void RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, size_t start_node,bool& stop);
        bool ProcessNeighbors(std::unordered_set<uint64_t>& neighbors, int& current_node, uint64_t start_node, std::vector<uint64_t>& path, std::map<uint64_t, int>& lock, std::vector<std::unordered_set<uint64_t>>& stack, std::vector<int>& backtrack_lengths, std::vector<std::vector<uint64_t>>& cycles, int& counter,bool& stop);
        void Backtrack(std::vector<std::unordered_set<uint64_t>>& stack, std::vector<uint64_t>& path, std::vector<int>& backtrack_lengths, std::map<uint64_t, int>& lock, uint64_t start_node,bool& stop);
        vector<vector<uint64_t>> FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<unordered_set<uint64_t>> stack, vector<int> backtrack_lengths,bool& stop);
        vector<vector<uint64_t>> FindCycleUtil(uint64_t startnode);
        bool DepthLevelSearch(uint64_t start, uint64_t target, int limit, int& reached_depth);
        size_t ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked);
        int FindCycles();
        ~CycleFinder();
};  
#endif