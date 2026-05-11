#ifndef PIPELINE_H
#define PIPELINE_H

#include "settings.h"
#include "sdbg/sdbg.h"
#include "sdbg_build.h"
#include "cycle_finder.h"
#include "post_processing.h"
#include <iostream>
#include <filesystem>

inline void step_build_graph(Settings& settings) {
    SDBGBuild sdbg_build(settings);
}

inline void step_load_graph(Settings& settings, SDBG& sdbg) {
    settings.graph_folder += "/graph";
    std::cout << "Graph folder: " << settings.graph_folder << std::endl;
    sdbg.LoadFromFile(settings.graph_folder.c_str());
    std::cout << "Loaded the graph" << std::endl;
    settings.sdbg = &sdbg;
}

inline void step_find_cycles(Settings& settings) {
    std::cout << "FBCE START:" << std::endl;
    CycleFinder cycle_finder(settings);
    std::cout << "Number of nodes in results: " << cycle_finder.results.size() << std::endl;
}

inline void step_post_process(Settings& settings) {
    std::cout << "POST PROCESSING START:" << std::endl;
    PostProcessor processor(settings);
    processor.run_analysis();
    std::cout << "Saved in: " << settings.output_folder << std::endl;
}

inline void step_cleanup(const Settings& settings) {
    namespace fs = std::filesystem;
    // Remove the SDBG graph folder (large binary files, no longer needed)
    fs::path graph_dir = fs::path(settings.output_folder) / "graph";
    std::error_code ec;
    fs::remove_all(graph_dir, ec);
    if (ec) std::cerr << "Warning: could not remove graph folder: " << ec.message() << std::endl;
    // Remove the intermediate cycles file
    fs::path cycles_file = fs::path(settings.output_folder) / "cycles.txt";
    fs::remove(cycles_file, ec);
    if (ec) std::cerr << "Warning: could not remove cycles.txt: " << ec.message() << std::endl;
    // Remove the empty cycles subfolder
    fs::path cycles_dir = fs::path(settings.output_folder) / "cycles";
    fs::remove(cycles_dir, ec); // only removes if empty
}

#endif // PIPELINE_H
