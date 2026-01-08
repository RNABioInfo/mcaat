#include <iostream>
#include <cstring>
#include "sdbg/sdbg.h"
#include "amino_acidator.h"

using namespace std;

int main(int argc, char** argv) {
    cout << "====================================" << endl;
    cout << "Testing AminoAcidator" << endl;
    cout << "====================================" << endl;
    
    // Load the graph
    SDBG sdbg;
    string graph_path = "/home/alex/git/mcaat/build/mcaat_run_2025-11-22_00-05-44/graph/graph";
    
    if (argc > 1) {
        graph_path = argv[1];
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
    
    // Create AminoAcidator instance
    AminoAcidator acidator(sdbg);
    
    // Find a valid starting node with diverse outgoing edges
    uint64_t start_node = 1000000;
    bool found_valid_start = false;
    
    for (uint64_t i = 100000; i < min(sdbg.size(), (uint64_t)200000); ++i) {
        if (sdbg.IsValidEdge(i) && sdbg.EdgeOutdegree(i) >= 2) {
            start_node = i;
            found_valid_start = true;
            break;
        }
    }
    
    if (!found_valid_start) {
        cout << "Could not find a valid start node!" << endl;
        return 1;
    }
    
    cout << "Starting from node: " << start_node << endl;
    cout << "Node outdegree: " << sdbg.EdgeOutdegree(start_node) << endl;
    
    // Debug: Check the actual sequence
    vector<uint8_t> test_seq(sdbg.k());
    sdbg.GetLabel(start_node, test_seq.data());
    cout << "Start node sequence (raw): ";
    for (uint32_t i = 0; i < sdbg.k(); ++i) {
        cout << (int)test_seq[i] << " ";
    }
    cout << endl;
    
    cout << "Start node sequence (DNA): ";
    for (uint32_t i = 0; i < sdbg.k(); ++i) {
        switch (test_seq[i]) {
            case 1: cout << 'A'; break;
            case 2: cout << 'C'; break;
            case 3: cout << 'G'; break;
            case 4: cout << 'T'; break;
            default: cout << 'N'; break;
        }
    }
    cout << endl;
    
    // Test outgoing edges
    uint64_t test_outgoings[4];
    int test_outdegree = sdbg.OutgoingEdges(start_node, test_outgoings);
    cout << "Outgoing edges: ";
    for (int i = 0; i < test_outdegree; ++i) {
        cout << test_outgoings[i] << " ";
    }
    cout << endl << endl;
    
    // Run beam search with different parameters
    int beam_width = 5;
    int max_depth = 50;
    
    cout << "Running beam search with:" << endl;
    cout << "  Beam width: " << beam_width << endl;
    cout << "  Max depth: " << max_depth << endl;
    cout << endl;
    
    auto start_time = chrono::high_resolution_clock::now();
    vector<AminoAcidPathInfo> paths = acidator.BeamSearchAminoAcids(start_node, beam_width, max_depth);
    auto end_time = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
    cout << "Search completed in " << duration.count() << " ms" << endl;
    cout << "Found " << paths.size() << " paths" << endl;
    cout << endl;
    
    // Display results for top paths
    int num_to_display = min(3, (int)paths.size());
    cout << "Displaying top " << num_to_display << " paths:" << endl;
    cout << "-----------------------------------" << endl;
    
    for (int i = 0; i < num_to_display; ++i) {
        const auto& path = paths[i];
        cout << "\nPath " << (i + 1) << ":" << endl;
        cout << "  Nodes in path: " << path.node_path.size() << endl;
        cout << "  Amino acids: " << path.amino_acids.size() << endl;
        cout << "  Total score: " << path.total_score << endl;
        
        // Display amino acid sequence
        cout << "  AA sequence: ";
        for (const auto& aa : path.amino_acids) {
            cout << aa;
        }
        cout << endl;
        
        // Reconstruct the FULL DNA sequence to verify
        string full_dna_seq;
        for (size_t j = 0; j < path.node_path.size(); ++j) {
            vector<uint8_t> node_seq(sdbg.k());
            sdbg.GetLabel(path.node_path[j], node_seq.data());
            if (j == 0) {
                // First node: add all bases
                for (uint32_t idx = 0; idx < sdbg.k(); ++idx) {
                    switch (node_seq[idx]) {
                        case 1: full_dna_seq += 'A'; break;
                        case 2: full_dna_seq += 'C'; break;
                        case 3: full_dna_seq += 'G'; break;
                        case 4: full_dna_seq += 'T'; break;
                        default: full_dna_seq += 'N'; break;
                    }
                }
            } else {
                // Subsequent nodes: only last base (k-mer overlap)
                switch (node_seq[sdbg.k() - 1]) {
                    case 1: full_dna_seq += 'A'; break;
                    case 2: full_dna_seq += 'C'; break;
                    case 3: full_dna_seq += 'G'; break;
                    case 4: full_dna_seq += 'T'; break;
                    default: full_dna_seq += 'N'; break;
                }
            }
        }
        cout << "  Full DNA: " << full_dna_seq.substr(0, min((size_t)100, full_dna_seq.length())) << endl;
        cout << "  Accum DNA: " << path.dna_sequence.substr(0, min((size_t)100, path.dna_sequence.length())) << endl;
        
        // Verify they match
        if (full_dna_seq != path.dna_sequence) {
            cout << "  WARNING: DNA sequences don't match!" << endl;
            cout << "  Full length: " << full_dna_seq.length() << " vs Accum length: " << path.dna_sequence.length() << endl;
        }
        
        // Manually verify codon translation
        cout << "  Last 12bp DNA: " << path.dna_sequence.substr(path.dna_sequence.length() - 12) << endl;
        cout << "  Last 4 codons: ";
        for (int c = 4; c >= 1; c--) {
            size_t start_pos = path.dna_sequence.length() - (c * 3);
            if (start_pos < path.dna_sequence.length()) {
                string codon = path.dna_sequence.substr(start_pos, 3);
                cout << codon << "  ";
            }
        }
        cout << endl;
        
        // Display first few node IDs
        cout << "  Node path (first 10): ";
        int nodes_to_show = min(10, (int)path.node_path.size());
        for (int j = 0; j < nodes_to_show; ++j) {
            cout << path.node_path[j];
            if (j < nodes_to_show - 1) cout << " -> ";
        }
        if (path.node_path.size() > 10) {
            cout << " ... (" << (path.node_path.size() - 10) << " more)";
        }
        cout << endl;
        
        // Display first few scores
        cout << "  Scores (first 5): ";
        int scores_to_show = min(5, (int)path.scores.size());
        for (int j = 0; j < scores_to_show; ++j) {
            cout << path.scores[j];
            if (j < scores_to_show - 1) cout << ", ";
        }
        if (path.scores.size() > 5) {
            cout << " ... (" << (path.scores.size() - 5) << " more)";
        }
        cout << endl;
    }
    
    cout << endl;
    cout << "====================================" << endl;
    cout << "Test completed successfully!" << endl;
    cout << "====================================" << endl;
    
    return 0;
}
