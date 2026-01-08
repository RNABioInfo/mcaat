#include <iostream>
#include <filesystem>
#include "sdbg_build.h"
#include "settings.h"

int main() {
    Settings settings;
    
    // Set up basic parameters for graph building
    settings.input_files = "/home/alex/git/mcaat/build/sim_reads.fastq";
    settings.output_folder = "sim_graph_output";
    settings.graph_folder = "sim_graph_output/graph";
    settings.cycles_folder = "sim_graph_output/cycles";
    settings.threads = 13;
    settings.ram = 16.0; // 16GB
    
    // Create directories
    try {
        std::filesystem::create_directories(settings.output_folder);
        std::filesystem::create_directories(settings.graph_folder);
        std::filesystem::create_directories(settings.cycles_folder);
        std::cout << "Created output directories" << std::endl;
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating directories: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "Building graph from: " << settings.input_files << std::endl;
    std::cout << "Output to: " << settings.graph_folder << std::endl;
    std::cout << "Threads: " << settings.threads << std::endl;
    
    // Build the graph
    SDBGBuild sdbg_build(settings);
    
    std::cout << "Graph built successfully!" << std::endl;
    std::cout << "Graph location: " << settings.graph_folder << "/graph" << std::endl;
    
    return 0;
}
