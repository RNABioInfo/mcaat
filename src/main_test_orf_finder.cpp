/**
 * @file main_test_orf_finder.cpp
 * @brief Test program for ORF finder with SDBG
 * 
 * Builds a de Bruijn graph from paired-end FASTQ files and tests the ORF finder
 * to locate start/stop codons at specific distances from a known repeat k-mer.
 * 
 * Compile (macOS - without OpenMP):
 * g++ -std=c++17 -O3 -march=native \
 *     -Xpreprocessor -fopenmp \
 *     -I./include -I./libs/megahit/src -I./libs/kseqpp/include \
 *     -I/usr/local/opt/libomp/include \
 *     -L/usr/local/opt/libomp/lib \
 *     -o test_orf_finder \
 *     src/main_test_orf_finder.cpp src/orf_finder.cpp src/sdbg_build.cpp \
 *     libs/megahit/src/sdbg/sdbg_meta.cpp \
 *     libs/megahit/src/sdbg/sdbg_raw_content.cpp \
 *     libs/megahit/src/sdbg/sdbg_writer.cpp \
 *     libs/megahit/src/sorting/kmer_counter.cpp \
 *     libs/megahit/src/sorting/read_to_sdbg_s1.cpp \
 *     libs/megahit/src/sorting/read_to_sdbg_s2.cpp \
 *     libs/megahit/src/sorting/seq_to_sdbg.cpp \
 *     libs/megahit/src/utils/options_description.cpp \
 *     libs/megahit/src/sorting/base_engine.cpp \
 *     libs/megahit/src/sorting/kmsort_selector.cpp \
 *     libs/megahit/src/sequence/io/fastx_reader.cpp \
 *     libs/megahit/src/sequence/io/sequence_lib.cpp \
 *     libs/megahit/src/sequence/io/paired_fastx_reader.cpp \
 *     -lomp -lz
 * 
 * Usage:
 * ./test_orf_finder <reads_R1.fastq> <reads_R2.fastq> <output_dir>
 */

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <cstring>
#include "sdbg_build.h"
#include "settings.h"
#include "orf_finder.h"
#include "sdbg/sdbg.h"

namespace fs = std::filesystem;

// Target k-mer to find (repeat sequence)
const std::string TARGET_KMER = "GTTTCCATTAATTCCACTTCTTA";

// Helper: Convert k-mer string to node ID, handling k-mer length mismatches
uint64_t kmer_to_node_id(const SDBG& sdbg, const std::string& kmer) {
    const uint32_t k = sdbg.k();
    
    // If k-mer is longer than graph k, try all sub-kmers
    if (kmer.size() > k) {
        std::cout << "  K-mer (" << kmer.size() << " bp) is longer than graph k=" << k << std::endl;
        std::cout << "  Trying all " << k << "-mers within the sequence..." << std::endl;
        
        for (size_t start = 0; start <= kmer.size() - k; ++start) {
            std::string sub_kmer = kmer.substr(start, k);
            std::vector<uint8_t> seq(k);
            bool valid = true;
            
            for (uint32_t i = 0; i < k; ++i) {
                char c = sub_kmer[i];
                if (c == 'A') seq[i] = 1;
                else if (c == 'C') seq[i] = 2;
                else if (c == 'G') seq[i] = 3;
                else if (c == 'T') seq[i] = 4;
                else {
                    valid = false;
                    break;
                }
            }
            
            if (valid) {
                int64_t result = sdbg.IndexBinarySearch(seq.data());
                if (result != -1) {
                    std::cout << "  ✓ Found sub-kmer at position " << start << ": " << sub_kmer << std::endl;
                    return static_cast<uint64_t>(result);
                }
            }
        }
        
        std::cout << "  ✗ No matching sub-kmer found in graph" << std::endl;
        return SDBG::kNullID;
    }
    
    // Exact k-mer match
    if (kmer.size() != k) {
        std::cerr << "  Warning: k-mer size " << kmer.size() 
                  << " does not match k=" << k << std::endl;
        return SDBG::kNullID;
    }
    
    std::vector<uint8_t> seq(k);
    for (uint32_t i = 0; i < k; ++i) {
        char c = kmer[i];
        if (c == 'A') seq[i] = 1;
        else if (c == 'C') seq[i] = 2;
        else if (c == 'G') seq[i] = 3;
        else if (c == 'T') seq[i] = 4;
        else {
            std::cerr << "  Invalid nucleotide: " << c << std::endl;
            return SDBG::kNullID;
        }
    }
    
    int64_t result = sdbg.IndexBinarySearch(seq.data());
    return (result == -1) ? SDBG::kNullID : static_cast<uint64_t>(result);
}

// Helper: Get sequence from node
std::string node_to_sequence(const SDBG& sdbg, uint64_t node_id) {
    const uint32_t k = sdbg.k();
    std::vector<uint8_t> seq(k);
    sdbg.GetLabel(node_id, seq.data());
    
    std::string result;
    result.reserve(k);
    for (uint32_t i = 0; i < k; ++i) {
        switch(seq[i]) {
            case 1: result += 'A'; break;
            case 2: result += 'C'; break;
            case 3: result += 'G'; break;
            case 4: result += 'T'; break;
            default: result += 'N'; break;
        }
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <reads_R1.fastq> <reads_R2.fastq> <output_dir>" << std::endl;
        std::cerr << "\nExample:" << std::endl;
        std::cerr << "  " << argv[0] << " reads_R1.fastq reads_R2.fastq test_output" << std::endl;
        return 1;
    }
    
    std::string reads_r1 = argv[1];
    std::string reads_r2 = argv[2];
    std::string output_dir = argv[3];
    
    // Remove trailing slash
    while (!output_dir.empty() && output_dir.back() == '/') {
        output_dir.pop_back();
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << "ORF Finder Test - Paired-End Reads" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Input R1: " << reads_r1 << std::endl;
    std::cout << "Input R2: " << reads_r2 << std::endl;
    std::cout << "Output dir: " << output_dir << std::endl;
    std::cout << "Target k-mer: " << TARGET_KMER << " (" << TARGET_KMER.size() << " bp)" << std::endl;
    std::cout << std::endl;
    
    // Step 1: Build graph
    std::cout << "[1/4] Building de Bruijn graph from paired-end reads..." << std::endl;
    
    Settings settings;
    settings.input_files = reads_r1 + "," + reads_r2;  // Paired-end format
    settings.output_folder = output_dir;
    settings.graph_folder = output_dir + "/graph";
    settings.cycles_folder = output_dir + "/cycles";
    settings.threads = 4;
    settings.ram = 4.0;
    
    try {
        fs::create_directories(settings.output_folder);
        fs::create_directories(settings.graph_folder);
        fs::create_directories(settings.cycles_folder);
    } catch (const std::exception& e) {
        std::cerr << "Error creating directories: " << e.what() << std::endl;
        return 1;
    }
    
    SDBGBuild sdbg_build(settings);
    std::cout << "  Graph built at: " << settings.graph_folder << "/graph" << std::endl;
    std::cout << std::endl;
    
    // Step 2: Load graph
    std::cout << "[2/4] Loading graph..." << std::endl;
    
    SDBG sdbg;
    std::string graph_path = settings.graph_folder + "/graph";
    char* cstr = new char[graph_path.length() + 1];
    std::strcpy(cstr, graph_path.c_str());
    
    sdbg.LoadFromFile(cstr);
    delete[] cstr;
    
    std::cout << "  Graph size: " << sdbg.size() << " nodes" << std::endl;
    std::cout << "  K-mer size (k): " << sdbg.k() << std::endl;
    std::cout << std::endl;
    
    // Step 3: Find the target k-mer node
    std::cout << "[3/4] Searching for target k-mer in graph..." << std::endl;
    std::cout << "  Target: " << TARGET_KMER << std::endl;
    
    uint64_t repeat_node = kmer_to_node_id(sdbg, TARGET_KMER);
    
    if (repeat_node == SDBG::kNullID) {
        std::cerr << "\n✗ ERROR: Target k-mer not found in graph!" << std::endl;
        std::cerr << "  The k-mer '" << TARGET_KMER << "' does not exist in the assembled graph." << std::endl;
        std::cerr << "  This could mean:" << std::endl;
        std::cerr << "    - The sequence is not present in the input reads" << std::endl;
        std::cerr << "    - The k-mer was filtered during graph construction" << std::endl;
        std::cerr << "    - The k-mer length doesn't match the graph k-value" << std::endl;
        return 1;
    }
    
    std::string repeat_seq = node_to_sequence(sdbg, repeat_node);
    int mult = sdbg.EdgeMultiplicity(repeat_node);
    
    std::cout << "  ✓ Found target k-mer!" << std::endl;
    std::cout << "    Node ID: " << repeat_node << std::endl;
    std::cout << "    Sequence: " << repeat_seq << std::endl;
    std::cout << "    Multiplicity: " << mult << std::endl;
    std::cout << std::endl;
    
    // Step 4: Test ORF finder with different distance ranges
    std::cout << "[4/4] Running ORF finder on target k-mer..." << std::endl;
    std::cout << "========================================" << std::endl;
    
    ORFFinder orf_finder(sdbg);
    
    // Test multiple distance ranges for CRISPR-Cas genes
    std::vector<std::pair<int, int>> distance_ranges = {
        {50, 300},      // TINY genes (e.g., small Cas proteins)
        {100, 600},     // VERY-SMALL genes
        {200, 1000},    // SMALL genes
        {300, 1500},    // MEDIUM-SMALL genes
        {500, 2500},    // MEDIUM genes
    };
    
    bool found_any_orf = false;
    
    for (const auto& [min_dist, max_dist] : distance_ranges) {
        std::cout << "\nSearching distance range: [" << min_dist << ", " << max_dist << "] bp" << std::endl;
        
        auto orf_result = orf_finder.FindFirstORF(repeat_node, min_dist, max_dist);
        
        if (orf_result.has_value()) {
            const auto& orf = orf_result.value();
            std::string start_seq = node_to_sequence(sdbg, orf.start_node);
            std::string end_seq = node_to_sequence(sdbg, orf.end_node);
            
            std::cout << "  ✓ ORF FOUND!" << std::endl;
            std::cout << "    Start node: " << orf.start_node << std::endl;
            std::cout << "    Start sequence: " << start_seq << std::endl;
            std::cout << "    End node: " << orf.end_node << std::endl;
            std::cout << "    End sequence: " << end_seq << std::endl;
            std::cout << "    Distance from repeat: " << orf.distance << " bp" << std::endl;
            std::cout << "    Protein length: ~" << (orf.distance / 3) << " amino acids" << std::endl;
            
            found_any_orf = true;
        } else {
            std::cout << "  ✗ No ORF found in this range" << std::endl;
        }
    }
    
    // Summary
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Target k-mer: " << TARGET_KMER << std::endl;
    std::cout << "Node ID: " << repeat_node << std::endl;
    std::cout << "Multiplicity: " << mult << std::endl;
    std::cout << std::endl;
    
    if (found_any_orf) {
        std::cout << "✓ SUCCESS: ORF finder found potential coding sequences!" << std::endl;
        std::cout << "  The ORF finder successfully identified start and stop codons" << std::endl;
        std::cout << "  at various distances from the target repeat k-mer." << std::endl;
    } else {
        std::cout << "⚠ NO ORFs FOUND in any tested distance range." << std::endl;
        std::cout << "  This could mean:" << std::endl;
        std::cout << "    - No start codons (ATG/GTG/TTG) exist within the search ranges" << std::endl;
        std::cout << "    - No stop codons (TAA/TAG/TGA) follow the start codons" << std::endl;
        std::cout << "    - The graph connectivity prevents BFS traversal" << std::endl;
        std::cout << "  Try adjusting distance ranges or using different input data." << std::endl;
    }
    
    return 0;
}
