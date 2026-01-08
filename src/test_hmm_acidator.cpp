#include <iostream>
#include <cstring>
#include "sdbg/sdbg.h"
#include "amino_acidator.h"
#include "profile.h"

using namespace std;

int main(int argc, char** argv) {
    cout << "====================================" << endl;
    cout << "Testing AminoAcidator with HMM Profile" << endl;
    cout << "====================================" << endl;
    
    // Load the graph
    SDBG sdbg;
    string graph_path = "/home/alex/git/mcaat/build/mcaat_run_2025-11-22_00-05-44/graph/graph";
    string hmm_path = "hmm_test.hmm";
    uint64_t start_node = 0;
    bool use_custom_start = false;
    
    if (argc > 1) {
        graph_path = argv[1];
    }
    if (argc > 2) {
        hmm_path = argv[2];
    }
    if (argc > 3) {
        start_node = stoull(argv[3]);
        use_custom_start = true;
    }
    
    cout << "Loading graph from: " << graph_path << endl;
    char* cstr = new char[graph_path.length() + 1];
    strcpy(cstr, graph_path.c_str());
    
    sdbg.LoadFromFile(cstr);
    delete[] cstr;
    
    cout << "Graph loaded successfully!" << endl;
    cout << "Graph size: " << sdbg.size() << " nodes" << endl;
    cout << "K-mer size: " << sdbg.k() << endl;
    cout << endl;
    
    // Load HMM profile
    cout << "Loading HMM profile from: " << hmm_path << endl;
    Profile profile;
    if (!profile.LoadFromFile(hmm_path)) {
        cerr << "Failed to load profile!" << endl;
        return 1;
    }
    cout << endl;
    
    // Create AminoAcidator with HMM profile
    AminoAcidator acidator(sdbg, &profile);
    
    // Find a valid starting node
    bool found_valid_start = false;
    
    if (use_custom_start) {
        if (start_node < sdbg.size() && sdbg.IsValidEdge(start_node)) {
            found_valid_start = true;
            cout << "Using provided start node: " << start_node << endl;
        } else {
            cout << "Provided node " << start_node << " is not valid!" << endl;
            if (start_node >= sdbg.size()) {
                cout << "  Node ID exceeds graph size (" << sdbg.size() << ")" << endl;
            }
        }
    } else {
        start_node = 100012;
        for (uint64_t i = 100000; i < min(sdbg.size(), (uint64_t)200000); ++i) {
            if (sdbg.IsValidEdge(i) && sdbg.EdgeOutdegree(i) >= 2) {
                start_node = i;
                found_valid_start = true;
                break;
            }
        }
    }
    
    if (!found_valid_start) {
        cout << "Could not find a valid start node!" << endl;
        return 1;
    }
    
    cout << "Starting from node: " << start_node << endl;
    cout << "Node outdegree: " << sdbg.EdgeOutdegree(start_node) << endl;
    cout << endl;
    
    // Run beam search with HMM scoring
    // Need ~3 nodes per amino acid, HMM has 328 positions, so need ~1000 depth
    int beam_width = 10;
    int max_depth = 1200;  // Allow enough depth to reach full HMM length
    
    cout << "Running HMM-scored beam search with:" << endl;
    cout << "  Beam width: " << beam_width << endl;
    cout << "  Max depth: " << max_depth << endl;
    cout << "  HMM length: " << profile.GetLength() << " positions" << endl;
    cout << endl;
    
    auto start_time = chrono::high_resolution_clock::now();
    vector<AminoAcidPathInfo> paths = acidator.BeamSearchAminoAcids(start_node, beam_width, max_depth);
    auto end_time = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
    cout << "Search completed in " << duration.count() << " ms" << endl;
    cout << "Found " << paths.size() << " paths" << endl;
    cout << endl;
    
    // Display results for top paths
    int num_to_display = min(5, (int)paths.size());
    cout << "Displaying top " << num_to_display << " paths (sorted by HMM score):" << endl;
    cout << "-----------------------------------" << endl;
    
    for (int i = 0; i < num_to_display; ++i) {
        const auto& path = paths[i];
        cout << "\nPath " << (i + 1) << ":" << endl;
        cout << "  HMM score: " << path.total_score << endl;
        cout << "  Nodes in path: " << path.node_path.size() << endl;
        cout << "  Amino acids: " << path.amino_acids.size() << endl;
        cout << "  HMM position reached: " << path.hmm_position << " / " << profile.GetLength() << endl;
        
        // Display amino acid sequence
        cout << "  AA sequence: ";
        for (const auto& aa : path.amino_acids) {
            cout << aa;
        }
        cout << endl;
        
        // Compare with HMM consensus
        cout << "  HMM consensus (first " << min((int)path.amino_acids.size(), 50) << "): ";
        for (int j = 1; j <= min((int)path.amino_acids.size(), min(50, profile.GetLength())); ++j) {
            cout << profile.GetState(j).consensus;
        }
        cout << endl;
        
        // Show match quality for first few positions
        if (path.amino_acids.size() > 0) {
            cout << "  First 5 AA match scores: ";
            for (size_t j = 0; j < min((size_t)5, path.amino_acids.size()); ++j) {
                if (j + 1 <= (size_t)profile.GetLength()) {
                    char aa = path.amino_acids[j][0];
                    double score = profile.GetMatchEmission(j + 1, aa);
                    cout << aa << ":" << score << " ";
                }
            }
            cout << endl;
        }
    }
    
    cout << "\n====================================" << endl;
    cout << "Test completed successfully!" << endl;
    cout << "====================================" << endl;
    
    return 0;
}
