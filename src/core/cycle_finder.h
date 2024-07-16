#ifndef CYCLE_FINDER_H
#define CYCLE_FINDER_H
#include <../src/sdbg/sdbg.h>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <stack>
//#include <omp.h>
#include "../helpers/node.h"
#include "../helpers/generic_functions.h"
#include "../helpers/path_writer.h"
#include "../../libs/progressbar/include/progressbar.hpp"
#include <list>
#include <thread>
#include <mutex>
#include <malloc.h>
using namespace std;

class CycleFinder {
    private:
        uint8_t maximal_length;
        uint8_t minimal_length;
        SDBG& sdbg;
        string genome_name;
        uint8_t cluster_bounds;
        vector<bool> visited;

        void _GetOutgoings(uint64_t node, set<uint64_t>& outgoings_set);
        void _GetIncomings(uint64_t node, set<uint64_t>& incomings_set);
        bool _MarkNode(uint64_t neighbor, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths);
        void _RelaxLock(int backtrack_length, uint64_t v, vector<uint64_t>& path, map<uint64_t, int>& lock, vector<set<uint64_t>>& stack, vector<int>& backtrack_lengths);
    public:
        CycleFinder(SDBG& sdbg, int length_bound, int minimal_length, string genome_name);
        size_t RemoveLinearPathsBacktrack(uint64_t node);
        vector<vector<uint64_t>> FindCycle(uint64_t start_node, vector<uint64_t> path, map<uint64_t, int> lock, vector<set<uint64_t>> stack, vector<int> backtrack_lengths);
        vector<vector<uint64_t>> FindCycleUtil(uint64_t startnode);
        void ChunkStartNodes(map<int, vector<uint64_t>, std::greater<int>>& start_nodes_chunked);
        int FindCycles();
        ~CycleFinder();
};  
#endif