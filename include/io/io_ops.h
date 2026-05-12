#ifndef IO_OPS_H
#define IO_OPS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
#include <sdbg/sdbg.h>

#ifdef DEVELOP
#include <nlohmann/json.hpp>
using json = nlohmann::json;
#endif

struct io_ops {
#ifdef DEVELOP
    static std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>> read_cycles(const std::string& file_path);
    static void write_cycles(const std::string& file_path, const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles);
    static void write_nodes_gfa(const std::string& file_path, const SDBG& sdbg);
#endif
};

#endif // IO_OPS_H