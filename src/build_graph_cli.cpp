#include <iostream>
#include <filesystem>
#include "sdbg_build.h"
#include "settings.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_fastq> <output_dir>" << std::endl;
        return 1;
    }
    
    std::string input_fastq = argv[1];
    std::string output_dir = argv[2];
    
    // Remove trailing slash if present
    while (!output_dir.empty() && output_dir.back() == '/') {
        output_dir.pop_back();
    }
    
    Settings settings;
    
    // Set up parameters for graph building
    settings.input_files = input_fastq;
    settings.output_folder = output_dir;
    settings.graph_folder = output_dir + "/graph";
    settings.cycles_folder = output_dir + "/cycles";
    settings.threads = 13;
    settings.ram = 16.0; // 16GB
    
    // Create directories
    try {
        std::filesystem::create_directories(settings.output_folder);
        std::filesystem::create_directories(settings.graph_folder);
        std::filesystem::create_directories(settings.cycles_folder);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating directories: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "Building graph from: " << settings.input_files << std::endl;
    std::cout << "Output to: " << settings.graph_folder << std::endl;
    
    // Build the graph
    SDBGBuild sdbg_build(settings);
    
    std::cout << "Graph built successfully at: " << settings.graph_folder << "/graph" << std::endl;
    
    return 0;
}
