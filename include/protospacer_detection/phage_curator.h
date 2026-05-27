#ifndef PHAGE_CURATOR_H
#define PHAGE_CURATOR_H

#include <vector>
#include <cstdint>
#include <sdbg/sdbg.h>
#include "crispr_array/graph_generic_func.h"
#include <string>
#include <stack>
#include <functional>
#include <map>
#include <set>
#include "spoa/spoa.hpp"


using namespace std;

struct BeamPathInfo {
    std::vector<uint64_t> path;
    double total_mult;
};

class PhageCurator {

private:
    // Variables
    SDBG& sdbg;
    std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>> grouped_paths;
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> cycles;
    std::set<uint64_t> cycle_nodes;
    std::map<uint64_t, double> avg_spacers;
    
    string _FetchNodeLastBase(size_t node);
    string _FetchFirstNode(size_t node);
    string _ReconstructPath(const std::vector<uint64_t>& path);
public:
    //ctor
    PhageCurator(SDBG& sdbg);
    PhageCurator(SDBG& sdbg, const std::map<uint64_t, std::map<uint64_t, std::vector<std::vector<uint64_t>>>>& grouped_paths, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles);
    std::vector<std::string> reconstructed_sequences;
    std::vector<std::vector<uint64_t>> DepthLimitedPaths(uint64_t start, int lower, int higher, std::function<void(const std::vector<uint64_t>&)> path_callback = nullptr);
    std::vector<std::vector<uint64_t>> DepthLimitedPathsAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden);
    std::vector<std::vector<uint64_t>> ExtendFromGroupedPaths(int min_depth, int max_depth);
    std::vector<std::vector<uint64_t>> ProcessPathsFromFile(const std::string& file_path, size_t min_group_size);
    void ReconstructPaths(std::vector<std::vector<uint64_t>> paths);
    void WriteSequencesToFasta(const std::string& filename);
    std::vector<BeamPathInfo> BeamSearchPaths(uint64_t start, int length, int beam_width);
    std::vector<uint64_t> FindBestPathBeamFromGroupedPaths(int min_length, int beam_width);
    bool RevalidateAllNodesButSingleton();
    std::vector<std::vector<uint64_t>> BeamSearchPathsAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden, int beam_width, double min_mult, double max_mult, std::function<void(const std::vector<uint64_t>&)> path_callback = nullptr);
    std::vector<std::vector<uint64_t>> BeamSearchBackwardAvoiding(uint64_t start, int lower, int higher, const std::set<uint64_t>& forbidden, int beam_width, double min_mult, double max_mult);
    std::map<std::string,vector<string>> FindQualityPathsBeamSearchFromGroupedPaths(int min_length, int max_length, const std::string& filename, int beam_width);
    std::map<std::string,vector<string>> FindContiguousPathsFromGroupedPaths(int max_total_bp, const std::string& filename, int beam_width);
    std::string ComputeConsensusForCurrentGroup(vector<string> sequences);
    std::vector<vector<uint64_t>> GetTopPathsFromBeamPaths(const std::vector<std::vector<uint64_t>>& beam_paths,int max,int min,size_t top_n);

};

#endif // PHAGE_CURATOR_H