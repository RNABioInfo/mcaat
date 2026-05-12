#ifndef PIPELINE_H
#define PIPELINE_H

#include "core/settings.h"
#include "sdbg/sdbg.h"
#include "core/sdbg_build.h"
#include "crispr_array/cycle_finder.h"
#include "crispr_array/post_processing.h"
#include <iostream>
#include <filesystem>

inline void step_build_graph(Settings& settings) {
    SDBGBuild sdbg_build(settings);
}

inline void step_load_graph(Settings& settings, SDBG& sdbg) {
    // Use pre-built graph if provided via --graph, otherwise use the built graph subfolder
    std::string load_path = settings.graph_input.empty()
        ? settings.graph_folder + "/graph"
        : settings.graph_input;
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

inline void step_cleanup(const Settings& settings) {
    namespace fs = std::filesystem;
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
