#include <iostream>
#include <cstring>
#include "sdbg/sdbg.h"
#include "cas_gene_detector.h"
#include "profile.h"

using namespace std;

int main(int argc, char** argv) {
    cout << "====================================" << endl;
    cout << "Testing CasGeneDetector with HMM Profile" << endl;
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
    
    // Create CasGeneDetector with HMM profile
    CasGeneDetector detector(sdbg, &profile);
    
    // Validate starting node
    if (!use_custom_start) {
        cerr << "ERROR: No start node provided!" << endl;
        cerr << "Usage: " << argv[0] << " <graph_path> <hmm_path> <start_node>" << endl;
        return 1;
    }
    
    if (start_node >= sdbg.size() || !sdbg.IsValidEdge(start_node)) {
        cerr << "ERROR: Invalid start node " << start_node << endl;
        if (start_node >= sdbg.size()) {
            cerr << "  Node ID " << start_node << " exceeds graph size (" << sdbg.size() << ")" << endl;
        } else {
            cerr << "  Node is not a valid edge in the graph" << endl;
        }
        return 1;
    }
    
    cout << "Starting from node: " << start_node << endl;
    cout << "Node outdegree: " << sdbg.EdgeOutdegree(start_node) << endl;
    cout << endl;
    
    // Run beam search with HMM scoring
    // Each amino acid needs ~3 nucleotides (codon)
    // In de Bruijn graph, we traverse multiple nodes per position
    int beam_width = 10;
    int max_depth = profile.GetLength() * 3;  // HMM length * 3 bases per AA
    
    cout << "Running HMM-scored beam search with:" << endl;
    cout << "  Beam width: " << beam_width << endl;
    cout << "  Max depth: " << max_depth << " (HMM length " << profile.GetLength() << " * 3)" << endl;
    cout << "  HMM length: " << profile.GetLength() << " positions" << endl;
    cout << endl;
    
    auto start_time = chrono::high_resolution_clock::now();
    vector<AminoAcidPathInfo> paths = detector.BeamSearchAminoAcids(start_node, beam_width, max_depth);
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
        
        // Beam search already used Viterbi scoring
        double beam_viterbi_score = path.total_score;
        
        // Re-run Viterbi to show alignment details
        auto [viterbi_bit_score, alignment_path, hmm_end_pos] = profile.ViterbiAlign(path.amino_acids);
        int viterbi_matches = 0;
        for (char c : alignment_path) {
            if (c == 'M') viterbi_matches++;
        }
        
        cout << "\nPath " << (i + 1) << ":" << endl;
        cout << "  Viterbi score: " << viterbi_bit_score << " bits" << endl;
        cout << "  (Beam search ranked paths by this Viterbi score)" << endl;
        cout << "  Viterbi alignment: " << alignment_path.substr(0, min((size_t)50, alignment_path.size())) 
             << (alignment_path.size() > 50 ? "..." : "") << endl;
        cout << "  Nodes in path: " << path.node_path.size() << endl;
        cout << "  Amino acids: " << path.amino_acids.size() << endl;
        cout << "  HMM matches (Viterbi): " << viterbi_matches << " / " << profile.GetLength() << endl;
        
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
