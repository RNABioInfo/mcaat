#include "io/io_ops.h"
#ifdef DE

std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> io_ops::read_cycles(const std::string& file_path) {
    std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> cycles;
    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return cycles;
    }
    json j;
    file >> j;
    for (auto& [key, value] : j.items()) {
        uint64_t cycle_id = std::stoull(key);
        std::vector<std::vector<uint64_t>> cycle_vec;
        for (const auto& arr : value) {
            std::vector<uint64_t> inner_vec = arr.get<std::vector<uint64_t>>();
            cycle_vec.push_back(inner_vec);
        }
        cycles[cycle_id] = cycle_vec;
    }
    file.close();
    return cycles;
}

void io_ops::write_cycles(const std::string& file_path, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles) {
    json j;
    for (const auto& [cycle_id, cycle_vec] : cycles) {
        j[std::to_string(cycle_id)] = cycle_vec;
    }
    std::ofstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    }
    file << j.dump(4); // Dump with indentation of 4 spaces
    file.close();
}


void io_ops::write_nodes_gfa(const std::string& file_path, const SDBG& sdbg) {
    std::ofstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "Error opening GFA file: " << file_path << std::endl;
        return;
    }
    // Write segments (nodes)
    for (uint64_t edge_id = 0; edge_id < sdbg.size(); ++edge_id) {
        if (!sdbg.IsValidEdge(edge_id)) continue;
        // Get label as sequence (optional, here just node id as string)
        file << "S\t" << edge_id << "\t*\n";
    }
    // Write links (edges)
    for (uint64_t edge_id = 0; edge_id < sdbg.size(); ++edge_id) {
        if (!sdbg.IsValidEdge(edge_id)) continue;
        int outdeg = sdbg.EdgeOutdegree(edge_id);
        if (outdeg > 0) {
            std::vector<uint64_t> outgoings(outdeg);
            sdbg.OutgoingEdges(edge_id, outgoings.data());
            for (int i = 0; i < outdeg; ++i) {
                // Overlap is unknown, use 0M (no overlap info)
                file << "L\t" << edge_id << "\t+\t" << outgoings[i] << "\t+\t0M\n";
            }
        }
    }
    file.close();
}

#endif