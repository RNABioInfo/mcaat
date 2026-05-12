#ifndef PIPELINE_H
#define PIPELINE_H

#include "core/settings.h"
#include "sdbg/sdbg.h"
#include "core/sdbg_build.h"
#include "core/array_to_nodes.h"
#include "crispr_array/cycle_finder.h"
#include "crispr_array/post_processing.h"
#include "cas/crispr_postprocessor.h"
#include <iostream>
#include <filesystem>

inline void step_build_graph(Settings& settings) {
    SDBGBuild sdbg_build(settings);
}

inline void step_load_graph(Settings& settings, SDBG& sdbg) {
    namespace fs = std::filesystem;
    std::string load_path;
    if (settings.graph_input.empty()) {
        load_path = settings.graph_folder + "/graph";
    } else {
        // If the user passed a directory, auto-append the prefix
        load_path = fs::is_directory(settings.graph_input)
            ? settings.graph_input + "/graph"
            : settings.graph_input;
    }
    std::cout << "  ▸ Loading graph: " << load_path << std::endl;
    sdbg.LoadFromFile(load_path.c_str());
    std::cout << "  ▸ Graph loaded  (" << sdbg.size() << " nodes, k=" << sdbg.k() << ")" << std::endl;
    settings.sdbg = &sdbg;
}

inline void step_find_cycles(Settings& settings) {
    CycleFinder cycle_finder(settings);
    std::cout << "  ▸ Raw cycle entries: " << cycle_finder.results.size() << " start nodes" << std::endl;
}

inline void step_post_process(Settings& settings) {
    PostProcessor processor(settings);
    processor.run_analysis();
}

/**
 * @brief Map post-processed CRISPR arrays back to SDBG nodes.
 *
 * Reads all CRISPR_Arrays_*.txt files from settings.output_folder,
 * encodes each repeat/spacer sequence as k-mers, and looks up the
 * corresponding SDBG node IDs via IndexBinarySearch.  Returns a
 * vector of FilteredArray structs ready for CasWorkflow.
 */
inline std::vector<CRISPRPostProcessor::FilteredArray>
step_build_filtered_arrays(const Settings& settings, const SDBG& sdbg) {
    std::cout << "  ▸ Mapping CRISPR arrays to graph nodes..." << std::endl;
    auto arrays = BuildFilteredArraysFromDir(settings.output_folder, sdbg);
    std::cout << "  ▸ Mapped " << arrays.size() << " arrays" << std::endl;
    return arrays;
}

inline void step_cleanup(const Settings& settings) {    namespace fs = std::filesystem;
    fs::path graph_dir = fs::path(settings.output_folder) / "graph";
    std::error_code ec;
    fs::remove_all(graph_dir, ec);
    if (ec) std::cerr << "  Warning: could not remove graph folder: " << ec.message() << std::endl;
    fs::path cycles_file = fs::path(settings.output_folder) / "cycles.txt";
    fs::remove(cycles_file, ec);
    fs::path cycles_dir = fs::path(settings.output_folder) / "cycles";
    fs::remove(cycles_dir, ec);
    std::cout << "  ▸ Intermediate files removed" << std::endl;
}

#endif // PIPELINE_H
