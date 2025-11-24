#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cctype>
#include <unordered_map>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
#include <sys/sysinfo.h>
#elif defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#ifdef DEBUG
#include "io_ops.h"
#endif

#include "main_run_and_debug.h"
#include "cycle_finder.h"
#include "filters.h"
#include "post_processing.h"
#include "sdbg/sdbg.h"
#include "sdbg_build.h"
#include "settings.h"
#include "phage_curator.h"
#include "isolate_protospacers.h"

using namespace std;
namespace fs = std::filesystem;
#ifdef DEBUG
#pragma message("DEBUG is defined")
#else
#pragma message("DEBUG is NOT defined")
#endif
void print_usage(const char* program_name) {
    cout << "-------------------------------------------------------" << endl;
    cout << "\n";
    cout << "mCAAT - Metagenomic CRISPR Array Analysis Tool v. 0.4" << endl;
    cout << "\n";
    cout << "-------------------------------------------------------" << endl;
}

bool check_for_error(Settings& settings) {
    cout << "Step 1. Checking the inputs: " << endl;
    string erroneous_message = settings.print_settings();
    if (erroneous_message.empty()) {
        cout << "All inputs are correct. [✔]" << endl;
        return false;
    }
    cout << "Please check the following: " << erroneous_message << endl;
    return true;
}

// Get total system RAM in GB
double get_total_system_ram() {
    double total_bytes = 0.0;
    
#ifdef __linux__
    struct sysinfo mem_info;
    if (sysinfo(&mem_info) == 0) {
        total_bytes = static_cast<double>(mem_info.totalram) * mem_info.mem_unit;
    }
#elif defined(_WIN32)
    MEMORYSTATUSEX mem_info;
    mem_info.dwLength = sizeof(MEMORYSTATUSEX);
    if (GlobalMemoryStatusEx(&mem_info)) {
        total_bytes = static_cast<double>(mem_info.ullTotalPhys);
    }
#elif defined(__APPLE__)
    int mib[2] = {CTL_HW, HW_MEMSIZE};
    int64_t mem_size;
    size_t len = sizeof(mem_size);
    if (sysctl(mib, 2, &mem_size, &len, nullptr, 0) == 0) {
        total_bytes = static_cast<double>(mem_size);
    }
#endif
    
    // Convert bytes to GB
    return total_bytes / (1024.0 * 1024.0 * 1024.0);
}

Settings parse_arguments(int argc, char* argv[]) {
    vector<string> input_files_default;
    Settings settings;

    // Get timestamp for fallback folder naming
    string timestamp = settings.get_timestamp();

    // Pre-scan for --settings so file values can be used as defaults and CLI overrides them
    for (int j = 1; j < argc; ++j) {
        if (string(argv[j]) == "--settings" && j + 1 < argc) {
            string settings_file = argv[j + 1];
            if (!settings.LoadFromFile(settings_file))
                throw runtime_error("Error: could not load settings from " + settings_file);
            break;
        }
    }

     bool output_folder_provided = false;
 bool required_files_provided = false;
     bool input_files_from_settings = false;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];

        if (arg == "--help" || arg == "-h" || arg=="") {
            cout << "Usage: ./crispr_analyzer --input_files <file1> [file2] [options]\n"
                 << "\nRequired:\n"
                 << "  --input_files <file1> [file2]   One or two input FASTA/FASTQ files\n"
                 << "\nOptional:\n"
                 << "  --ram <amount>                  RAM to use (e.g., 4G, 500M). Default: 95% of system RAM\n"
                 << "  --threads <num>                 Number of threads. Default: CPU cores - 2\n"
                 << "  --output-folder <path>          Output directory. If not provided, a timestamped folder is created\n"
                 << "  --benchmark <file>              File containing expected crispr sequences line separated\n"
                 << "  --cycle-max-length <int>       Maximum cycle length to search (default in settings)\n"
                 << "  --cycle-min-length <int>       Minimum cycle length to search (default in settings)\n"
                 << "  --threshold-multiplicity <int> Minimum multiplicity threshold for start nodes (default in settings)\n"
                 << "  --low-abundance <true|false>   Enable low abundance mode for cycle filtering\n"
                 << "  --settings <path>              Path to a key=value settings file (overridden by CLI args)\n"
                 << "  --help, -h                      Show this help message\n";
            exit(0);
        }
       
        if (arg == "--input-files" || arg == "-i") {
            while (++i < argc && argv[i][0] != '-') {
                input_files_default.push_back(argv[i]);
                
            }
            --i;
            required_files_provided = true;
        } else if (arg == "--benchmark") {
            if (++i < argc) {
                settings.benchmark_file = argv[i];
            } else {
                throw runtime_error("Error: Missing value for --benchmark");
            }
            --i;
        } else if (arg == "--ram") {
            if (++i < argc) {
                string ram_input = argv[i];
                try {
                    double value = 0.0;
                    char unit = 'G';
                    size_t pos = ram_input.find_first_not_of("0123456789.");
                    if (pos != string::npos) {
                        value = stod(ram_input.substr(0, pos));
                        unit = toupper(ram_input[pos]);
                    } else {
                        value = stod(ram_input);
                    }

                    switch (unit) {
                        case 'B': settings.ram = value / (1024.0 * 1024.0 * 1024.0); break;
                        case 'K': settings.ram = value / (1024.0 * 1024.0); break;
                        case 'M': settings.ram = value / 1024.0; break;
                        case 'G': settings.ram = value; break;
                        default:
                            throw runtime_error("Error: Invalid RAM unit. Use B, K, M, or G.");
                    }

                    double total_ram = get_total_system_ram();
                    if (settings.ram < 1.0) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GB is too low (must be at least 1 GB)");
                    }
                    if (settings.ram > total_ram) {
                        throw runtime_error("Error: RAM value " + to_string(settings.ram) +
                                            " GB exceeds system total of " + to_string(total_ram) + " GB");
                    }
                } catch (const invalid_argument&) {
                    throw runtime_error("Error: Invalid RAM value provided: " + ram_input);
                }
            } else {
                throw runtime_error("Error: Missing value for --ram");
            }
        } else if (arg == "--threads") {
            if (++i < argc) {
                settings.threads = stoul(argv[i]);
                
            } else {
                throw runtime_error("Error: Missing value for --threads");
            }
        } else if (arg == "--output-folder" || arg == "--output_folder") {
            if (++i < argc) {
                settings.output_folder = string(argv[i]);
                output_folder_provided = true;
            } else {
                throw runtime_error("Error: Missing value for --output-folder");
            }
        }
        else if (arg == "--cycle-max-length") {
            if (++i < argc) {
                settings.cycle_finder_settings.cycle_max_length = stoi(argv[i]);
            } else {
                throw runtime_error("Error: Missing value for --cycle-max-length");
            }
        } else if (arg == "--cycle-min-length") {
            if (++i < argc) {
                settings.cycle_finder_settings.cycle_min_length = stoi(argv[i]);
            } else {
                throw runtime_error("Error: Missing value for --cycle-min-length");
            }
        } else if (arg == "--threshold-multiplicity") {
            if (++i < argc) {
                settings.cycle_finder_settings.threshold_multiplicity = stoull(argv[i]);
            } else {
                throw runtime_error("Error: Missing value for --threshold-multiplicity");
            }
        } else if (arg == "--low-abundance") {
            if (++i < argc) {
                string value = argv[i];
                std::transform(value.begin(), value.end(), value.begin(), ::tolower);
                settings.cycle_finder_settings.low_abundance = (value == "1" || value == "true" || value == "yes");
            } else {
                throw runtime_error("Error: Missing value for --low-abundance");
            }
        }
    }
    // If settings file provided input_files, allow it as equivalent to CLI --input-files
    if (input_files_default.empty() && !settings.input_files.empty()) {
        std::istringstream iss(settings.input_files);
        string token;
        while (iss >> token) {
            input_files_default.push_back(token);
        }
        required_files_provided = true;
        input_files_from_settings = true;
    }

    if (!required_files_provided && input_files_default.empty() && settings.input_files.empty()) {
        throw runtime_error("Error: No input files provided. Use --input-files <file1> [file2]");
    }
    // Set default output folder if not provided (only if not set by settings file)
    if (!output_folder_provided && settings.output_folder.empty()) {
        settings.output_folder = "mcaat_run_" + timestamp;
    }

    // Set subfolders and output file (but allow them to be set in settings file)
    if (settings.graph_folder.empty())
        settings.graph_folder = settings.output_folder + "/graph";
    if (settings.cycles_folder.empty())
        settings.cycles_folder = settings.output_folder + "/cycles";
    if (settings.output_file.empty())
        settings.output_file = settings.output_folder + "/CRISPR_Arrays.txt";

    // Debug output
    cout << "Output folder: " << settings.output_folder << endl;
        cout << "Graph folder: " << settings.graph_folder << endl;
        cout << "Cycles folder: " << settings.cycles_folder << endl;
        cout << "CycleFinder settings: max_length=" << settings.cycle_finder_settings.cycle_max_length
            << " min_length=" << settings.cycle_finder_settings.cycle_min_length
            << " threshold_mult=" << settings.cycle_finder_settings.threshold_multiplicity
            << " low_abundance=" << (settings.cycle_finder_settings.low_abundance ? "true" : "false")
            << " threads=" << settings.threads << endl;

    // Create directories
    try {
        fs::create_directories(settings.output_folder);
        fs::create_directories(settings.graph_folder);
        fs::create_directories(settings.cycles_folder);
        cout << "Created directories: " << settings.output_folder << ", "
             << settings.graph_folder << ", " << settings.cycles_folder << endl;
    } catch (const fs::filesystem_error& e) {
        throw runtime_error("Error: Could not create directories: " + string(e.what()));
    }

    // Validate input files
    if (input_files_default.size() < 1 || input_files_default.size() > 2) {
        throw runtime_error("Error: You must provide one or two input files.");
    }

    int count = 0;
    for (const auto& file : input_files_default) {
        if (!fs::exists(file)) {
            throw runtime_error("Error: Input file " + file + " does not exist.");
        }
        count++;
        // Only overwrite settings.input_files if the CLI actually provided files.
        // If input files were taken from the settings file, they are already in settings.input_files
        if (required_files_provided && !input_files_from_settings) {
            if (!settings.input_files.empty()) settings.input_files += " ";
            settings.input_files += file;
        }
    }

    if (settings.threads == 0) {
        settings.threads = thread::hardware_concurrency() - 2;
    }

    if (settings.ram == 0.0) {
        settings.ram = get_total_system_ram() * 0.95;
    }

    return settings;
}


string fetchNodeLabel(Settings& settings, uint64_t node) {
    std::string label;            
    uint8_t seq[settings.sdbg->k()];
    settings.sdbg->GetLabel(node, seq);
    for (int i = settings.sdbg->k() - 1; i >= 0; --i) label.append(1, "ACGT"[seq[i] - 1]);
    reverse(label.begin(), label.end());
    return label;
}

uint64_t findNodeFromKmer(Settings& settings, const std::string& kmer) {
    int k = settings.sdbg->k();
    if (static_cast<int>(kmer.size()) != k) {
        std::cerr << "Warning: kmer size " << kmer.size() << " does not match k=" << k << std::endl;
        return UINT64_MAX; // or handle error
    }
    std::vector<uint8_t> seq(k);
    for (int i = 0; i < k; ++i) {
        seq[i] = "ACGT"s.find(kmer[i]) + 1;
    }
    int64_t result = settings.sdbg->IndexBinarySearch(seq.data());
    return (result == -1) ? UINT64_MAX : static_cast<uint64_t>(result);
}

std::map<uint64_t, std::vector<std::vector<uint64_t>>> createRepeatToSpacerNodes(Settings& settings, const std::map<std::string, std::vector<std::string>>& systems_from_analyzer) {
    std::map<uint64_t, std::vector<std::vector<uint64_t>>> repeat_to_spacer_nodes;
    int k = settings.sdbg->k();
    for (const auto& [repeat, spacers] : systems_from_analyzer) {
        if (static_cast<int>(repeat.size()) < k) continue;
        std::string first_kmer = repeat.substr(0, k);
        uint64_t key_node = findNodeFromKmer(settings, first_kmer);
        if (key_node == UINT64_MAX) continue; // invalid
        std::vector<std::vector<uint64_t>> spacer_node_vectors;
        for (const auto& spacer : spacers) {
            std::vector<uint64_t> nodes;
            int L = spacer.size();
            for (int i = 0; i <= L - k; ++i) {
                std::string kmer = spacer.substr(i, k);
                uint64_t node = findNodeFromKmer(settings, kmer);
                if (node != UINT64_MAX) {
                    nodes.push_back(node);
                }
            }
            if (!nodes.empty()) {
                spacer_node_vectors.push_back(std::move(nodes));
            }
        }
        if (!spacer_node_vectors.empty()) {
            repeat_to_spacer_nodes[key_node] = std::move(spacer_node_vectors);
        }
    }
    return repeat_to_spacer_nodes;
}

#ifdef DEBUG
int main(int argc, char** argv) {
    // %% PARSE ARGUMENTS %%
    Settings settings = parse_arguments(argc, argv);
    string name_of_genome = "test";
    if (check_for_error(settings)){
        //tell the user which folder we are deleting
        cout<< "Folder " << settings.output_folder << " will be deleted due to errors." << endl;
        
        cout << "Do you want that folder to be removed? (y/n): ";
        char answer;
        cin >> answer;
        if (answer != 'y' && answer != 'Y') {
            cout << "Exiting the program." << endl;
            return 1;
        }
        cout << "Removing folder: " << settings.output_folder << endl;
        fs::remove_all(settings.output_folder); 
        return 1;
    }
    // %% PARSE ARGUMENTS %%

    // %% BUILD GRAPH %%
    SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
   
    // cycle finder max/min length are read from settings.cycle_finder_settings
    SDBG sdbg;
    vector<string> folders = {"/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-27_09-53-11/graph/graph","/vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-26-39/graph/graph",
    "/vol/d/development/git/mcaat_master/mcaat/build/mcaat_run_2025-10-24_12-27-47/graph/graph","/vol/d/data/real/graph"};
    string graph_folder_old = settings.graph_folder;///vol/d/development/git/mcaat_master/mcaat/_build/mcaat_run_2025-08-28_13-23-46
    settings.graph_folder=settings.graph_folder+"/graph";
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;
    settings.sdbg = &sdbg;

    // %% LOAD GRAPH %%
    
    delete[] cstr;
    string ending_seq = "CAGAGATAGAAATTATTTTTATTATACGTTTTTTTGT";
    vector<int64_t> end_nodes;
    //using findNodeFromKmer find all the node ids in the ending_seq
    for (size_t i = 0; i <= ending_seq.size() - settings.sdbg->k(); ++i) {
        string kmer = ending_seq.substr(i, sdbg.k());
        uint64_t node = findNodeFromKmer(settings, kmer);
        if (node != UINT64_MAX) {
            end_nodes.push_back(node);
        }
    }
    //print all the end nodes
    cout << "End nodes: ";
    for (const auto& node : end_nodes) {
        cout << node << " ";
    }
    cout << endl;
    // make some distribution of multiplicities of all nodes in the graph
    std::map<int, int> multiplicity_distribution;
    for (uint64_t node = 0; node < sdbg.size(); ++node) {
        int mult = sdbg.EdgeMultiplicity(node);
        multiplicity_distribution[mult]++;
    }
    cout << "Node Multiplicity Distribution:" << endl;
    //write distribution to file
    ofstream mult_file("node_multiplicities.txt");
    for (const auto& [mult, count] : multiplicity_distribution) {
        mult_file << "Multiplicity " << mult << ": " << count << " nodes" << endl;
    }

    mult_file.close();
    // for (uint64_t node = 0; node < sdbg.size(); ++node) {
    //     int mult = sdbg.EdgeMultiplicity(node);
    //     mult_file << node << "\t" << mult << "\n";
    // }
    // mult_file.close();
    // %% FBCE ALGORITHM %%
    cout << "FBCE FROM DEBUG START:" << endl;
    auto start_time = chrono::high_resolution_clock::now();
    CycleFinder cycle_finder(settings);
    // number_of_spacers_total not used; remove to avoid compiler warning
    auto cycles = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles.size() << endl;
    // %% FBCE ALGORITHM %%
    
    int number_of_spacers = 0;
    
    // %% FILTERS %%
    cout << "FILTERS START:" << endl;
    Filters filters(sdbg, cycles);
    auto  SYSTEMS = filters.ListArrays(number_of_spacers);
    cout<< "Number of spacers: " << number_of_spacers << " before cleaning"<<endl;
    // %% FILTERS %%

    //%% POST PROCESSING %%
    cout << "POST PROCESSING START:" << endl;
    CRISPRAnalyzer analyzer(SYSTEMS, settings.output_file);
    analyzer.run_analysis();
    cout << "Saved in: " << settings.output_file << endl;
    auto systems_from_analyzer = analyzer.getSystems();
    auto repeat_to_spacer_nodes = createRepeatToSpacerNodes(settings, systems_from_analyzer);
    cout << "Created repeat_to_spacer_nodes map with " << repeat_to_spacer_nodes.size() << " entries." << endl;
    //%% POST PROCESSING %%

    //%% PROTOSPACER ISOLATION %%
    IsolateProtospacers isolator(sdbg, repeat_to_spacer_nodes);
    pair<std::map<uint64_t,std::set<uint64_t>>,std::map<uint64_t,std::set<uint64_t>>> protospacer_nodes = isolator.getProtospacerNodes();
    auto grouped_paths_protospacers = isolator.DepthLimitedPathsFromInToOut(protospacer_nodes.first, protospacer_nodes.second, 50,1);
       isolator.WritePathsToFile(grouped_paths_protospacers, "grouped_paths_protospacers.txt");

    //auto grouped_paths_protospacers = isolator.ReadPathsFromFile("grouped_paths_protospacers.txt");
    
    //%% PROTOSPACER ISOLATION %%
    PhageCurator phage_curator(sdbg, grouped_paths_protospacers, cycles);

    // Find and output quality paths using DLS with multiplicity filters, writing live to file
    //phage_curator.FindQualityPathsDLSFromGroupedPaths(3000, 3010, "QualityPaths.fasta");  // min_length 3000, max_length arbitrary large
    vector<int> beam_widths = {50};
    // Find and output quality paths using beam search with multiplicity filters, writing live to file
    for (int beam_width : beam_widths) {
        string filename = "QualityPaths_BeamWidth" + to_string(beam_width) + ".fasta";
        phage_curator.FindQualityPathsBeamSearchFromGroupedPaths(3000, 3010, filename, beam_width);  // min_length 3000, max_length arbitrary large
    }
    //print multiplicities of nodes into a file

    //io_ops::write_nodes_gfa("output.gfa", sdbg);
    //io_ops::write_nodes_gfa("output.gfa", sdbg);
    

    //%% POST PROCESSING %%
    // %% DELETE THE GRAPH FOLDER %%
    //fs::remove_all(graph_folder_old);
    // %% DELETE THE GRAPH FOLDER %%          
    
}

#else
int main(int argc, char** argv) {
    // %% PARSE ARGUMENTS %%
    Settings settings = parse_arguments(argc, argv);
    if (check_for_error(settings)){
        cout << "Folder " << settings.output_folder;
        cout << " will be deleted due to errors." << endl;
        
        cout << "Do you want that folder to be removed? (y/n): ";
        char answer;
        cin >> answer;
        if (answer != 'y' && answer != 'Y') {
            cout << "Exiting the program." << endl;
            return 1;
        }
        cout << "Removing folder: " << settings.output_folder << endl;
        fs::remove_all(settings.output_folder); 
        return 1;
    }
    // %% PARSE ARGUMENTS %%

    // %% BUILD GRAPH %%
    SDBGBuild sdbg_build(settings);
    // %% BUILD GRAPH %%
    
    // %% LOAD GRAPH %%
    // cycle finder max/min length are read from settings.cycle_finder_settings
    SDBG sdbg;
    string graph_folder_old = settings.graph_folder;
    settings.graph_folder+="/graph";
    char * cstr = new char [settings.graph_folder.length()+1];
    std::strcpy (cstr, settings.graph_folder.c_str());
    cout << "Graph folder: " << cstr << endl;
    sdbg.LoadFromFile(cstr);
    cout << "Loaded the graph" << endl;
    settings.sdbg = &sdbg;
    // %% LOAD GRAPH %%

    // %% FBCE ALGORITHM %%
    cout << "FBCE START:" << endl;
    auto start_time = chrono::high_resolution_clock::now();
    CycleFinder cycle_finder(settings);
    // number_of_spacers_total unused
    auto cycles_map = cycle_finder.results;
    cout << "Number of nodes in results: " << cycles_map.size() << endl;
    // %% FBCE ALGORITHM %%
    
    auto cycles = cycles_map_to_cycles(cycles_map); // easier type to handle

    cout << "\n══════════════════════════════════════════════" << endl;
    cout << "🔸STEP 6: Finding relevant reads" << endl;
    cout << "══════════════════════════════════════════════" << endl;
    const auto reads = run_and_debug_finding_of_relevant_reads(
        cycles,
        settings,
        sdbg
    );

    cout << "\n══════════════════════════════════════════════" << endl;
    cout << "🔸STEP 7: Order the spacers" << endl;
    cout << "══════════════════════════════════════════════" << endl;
    const auto found_systems = run_and_debug_spacer_ordering(reads, sdbg, cycles);


    if (settings.benchmark_file != "") {
        cout << "\n══════════════════════════════════════════════" << endl;
        cout << "🔸STEP 8: Compare to ground of truth using benchmark file" << endl;
        cout << "══════════════════════════════════════════════" << endl;
        run_and_debug_benchmark_results(settings, found_systems);
    } else {
        cout << "\n══════════════════════════════════════════════" << endl;
        cout << "🔸STEP 8: Results" << endl;
        cout << "══════════════════════════════════════════════" << endl;
        run_and_debug_results(found_systems);
    }
    
    cout << "══════════════════════════════════════════════" << endl;

    // %% POST PROCESSING %%
    cout << "POST PROCESSING START:" << endl;
    unordered_map<string, vector<string>> all_systems;
    for (const auto& [_sequence, repeat, spacers, _conf_a, _conf_b] : found_systems) {
        all_systems[repeat] = spacers;
    }
    CRISPRAnalyzer analyzer(all_systems, settings.output_file);
    analyzer.run_analysis();
    cout << "Saved in: " << settings.output_file << endl;
    // %% POST PROCESSING %%

    // %% DELETE THE GRAPH FOLDER %%
    try {
        fs::remove_all(settings.graph_folder);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Warning: Could not remove graph folder: " << e.what() << std::endl;
    }
    // %% DELETE THE GRAPH FOLDER %%            
}
#endif