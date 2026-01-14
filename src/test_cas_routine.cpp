#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include "sdbg/sdbg.h"
#include "cas_workflow.h"
#include "orf_finder.h"

using namespace std;

/**
 * Find node ID for a given k-mer sequence in the SDBG
 */
uint64_t FindKmerNodeID(SDBG& sdbg, const string& kmer) {
    if (kmer.length() != (size_t)sdbg.k()) {
        cerr << "Error: K-mer length (" << kmer.length() 
             << ") doesn't match graph k (" << sdbg.k() << ")" << endl;
        return SDBG::kNullID;
    }

    // Convert DNA sequence to encoded format (A=1, C=2, G=3, T=4)
    uint8_t* encoded = new uint8_t[sdbg.k()];
    for (uint32_t i = 0; i < sdbg.k(); i++) {
        switch(kmer[i]) {
            case 'A': case 'a': encoded[i] = 1; break;
            case 'C': case 'c': encoded[i] = 2; break;
            case 'G': case 'g': encoded[i] = 3; break;
            case 'T': case 't': encoded[i] = 4; break;
            default: 
                cerr << "Invalid base: " << kmer[i] << endl;
                delete[] encoded;
                return SDBG::kNullID;
        }
    }
    
    uint64_t node_id = sdbg.IndexBinarySearch(encoded);
    delete[] encoded;
    
    return node_id;
}

/**
 * Get reverse complement of DNA sequence
 */
string ReverseComplement(const string& seq) {
    string rc = seq;
    for (size_t i = 0; i < rc.length(); i++) {
        switch(rc[i]) {
            case 'A': case 'a': rc[i] = 'T'; break;
            case 'T': case 't': rc[i] = 'A'; break;
            case 'C': case 'c': rc[i] = 'G'; break;
            case 'G': case 'g': rc[i] = 'C'; break;
        }
    }
    reverse(rc.begin(), rc.end());
    return rc;
}

int main(int argc, char** argv) {
    cout << "============================================" << endl;
    cout << "CAS WORKFLOW TEST - Operon Detection" << endl;
    cout << "============================================" << endl;
    cout << endl;
    
    // Parse command line arguments
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <graph_path> <profiles_dir> <repeat_sequence> [init_min] [init_max]" << endl;
        cerr << endl;
        cerr << "Arguments:" << endl;
        cerr << "  graph_path      - Path to SDBG graph file (e.g., graph/graph)" << endl;
        cerr << "  profiles_dir    - Directory containing HMM profile files" << endl;
        cerr << "  repeat_sequence - K-mer sequence from repeat region" << endl;
        cerr << "  init_min        - Optional: Min initial search distance (default: 0)" << endl;
        cerr << "  init_max        - Optional: Max initial search distance (default: 1000)" << endl;
        cerr << endl;
        return 1;
    }
    
    string graph_path = argv[1];
    string profiles_dir = argv[2];
    string repeat_kmer = argv[3];
    int init_min = (argc > 4) ? stoi(argv[4]) : 50;
    int init_max = (argc > 5) ? stoi(argv[5]) : 5000;
    
    cout << "Configuration:" << endl;
    cout << "  Graph path:      " << graph_path << endl;
    cout << "  Profiles dir:    " << profiles_dir << endl;
    cout << "  Repeat k-mer:    " << repeat_kmer << endl;
    cout << "  Initial search:  [" << init_min << ", " << init_max << "] bp" << endl;
    cout << endl;
    
    // Load SDBG graph
    cout << "Loading SDBG graph..." << endl;
    SDBG sdbg;
    char* graph_cstr = new char[graph_path.length() + 1];
    strcpy(graph_cstr, graph_path.c_str());
    sdbg.LoadFromFile(graph_cstr);
    delete[] graph_cstr;
    
    cout << "Graph loaded successfully!" << endl;
    cout << "  Graph size: " << sdbg.size() << " nodes" << endl;
    cout << "  K-mer size: " << sdbg.k() << endl;
    cout << endl;
    
    // Find repeat node ID from k-mer sequence
    cout << "Finding repeat node from k-mer sequence..." << endl;
    uint64_t repeat_node = FindKmerNodeID(sdbg, repeat_kmer);
    
    if (repeat_node == SDBG::kNullID) {
        cout << "K-mer not found in graph. Trying reverse complement..." << endl;
        string rc_kmer = ReverseComplement(repeat_kmer);
        cout << "  RC k-mer: " << rc_kmer << endl;
        repeat_node = FindKmerNodeID(sdbg, rc_kmer);
        
        if (repeat_node == SDBG::kNullID) {
            cerr << "ERROR: Repeat k-mer not found in graph (tried both strands)!" << endl;
            return 1;
        }
    }
    
    cout << "✓ Found repeat node: " << repeat_node << endl;
    
    // Verify the node
    if (!sdbg.IsValidEdge(repeat_node)) {
        cerr << "ERROR: Node " << repeat_node << " is not a valid edge!" << endl;
        return 1;
    }
    
    cout << "  Node outdegree: " << sdbg.EdgeOutdegree(repeat_node) << endl;
    cout << endl;
    
    // Create CAS workflow
    cout << "Initializing CAS workflow..." << endl;
    CasWorkflow workflow(sdbg, profiles_dir);
    
    // Set parameters
    workflow.SetParameters(
        init_min,        // initial_min (default 50bp)
        init_max,        // initial_max (default 5000bp)
        100,             // subsequent_search_distance (max 100bp between genes)
        41591,           // max_total_length
        100,             // beam_width (per spec)
        500              // max_search_depth (node depth fallback)
    );
    
    cout << "Workflow initialized with parameters:" << endl;
    cout << "  Subsequent search: 0-100 bp (handles 1-4bp overlaps and gaps)" << endl;
    cout << "  Max operon length: 41591 bp" << endl;
    cout << "  Beam width: 10" << endl;
    cout << "  Max search depth: 500" << endl;
    cout << "  Max nodes traversed: 5000" << endl;
    cout << endl;
    
    // Run CAS operon detection
    cout << "Starting CAS operon detection..." << endl;
    cout << "============================================" << endl;
    
    CasOperonResult result = workflow.DetectCasOperon(repeat_node);
    
    cout << "============================================" << endl;
    cout << endl;
    
    // Display results
    cout << "RESULTS:" << endl;
    cout << "--------" << endl;
    cout << "Genes detected: " << result.genes.size() << endl;
    cout << "Total operon length: " << result.total_length << " bp" << endl;
    cout << endl;
    
    if (result.genes.empty()) {
        cout << "No CAS genes detected." << endl;
        return 0;
    }
    
    cout << "Detected CAS Genes:" << endl;
    cout << "-------------------" << endl;
    
    for (size_t i = 0; i < result.genes.size(); i++) {
        const auto& gene = result.genes[i];
        
        cout << "Gene #" << (i + 1) << ":" << endl;
        cout << "  Name:             " << gene.gene_name << endl;
        cout << "  Start node:       " << gene.start_node << endl;
        cout << "  End node:         " << gene.end_node << endl;
        cout << "  Distance from R:  " << gene.distance_from_repeat << " bp" << endl;
        cout << "  Gene length:      " << gene.gene_length << " bp" << endl;
        cout << "  Score:            " << gene.score << endl;
        cout << "  Amino acids:      " << gene.amino_acids.size() << " AA" << endl;
        
        // Show first few amino acids
        if (!gene.amino_acids.empty()) {
            cout << "  AA sequence:      ";
            int num_to_show = min(20, (int)gene.amino_acids.size());
            for (int j = 0; j < num_to_show; j++) {
                cout << gene.amino_acids[j];
            }
            if (gene.amino_acids.size() > 20) {
                cout << "... (" << gene.amino_acids.size() << " total)";
            }
            cout << endl;
        }
        
        cout << "  Node path length: " << gene.orf_node_path.size() << " nodes" << endl;
        cout << endl;
    }
    
    cout << "Test completed successfully!" << endl;
    
    return 0;
}
