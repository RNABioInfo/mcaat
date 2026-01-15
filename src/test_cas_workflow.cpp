#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include <iomanip>
#include "sdbg/sdbg.h"
#include "cas_workflow.h"

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
    cout << "   CAS WORKFLOW TEST - Gene Detection" << endl;
    cout << "============================================" << endl;
    cout << endl;
    
    // Parse command line arguments
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <graph_path> <profiles_dir> <repeat_kmer>" << endl;
        cerr << endl;
        cerr << "Arguments:" << endl;
        cerr << "  graph_path    - Path to SDBG graph file (e.g., output/graph/graph)" << endl;
        cerr << "  profiles_dir  - Directory containing HMM profile files (e.g., profiles/)" << endl;
        cerr << "  repeat_kmer   - K-mer sequence from CRISPR repeat region" << endl;
        cerr << endl;
        cerr << "Example:" << endl;
        cerr << "  " << argv[0] << " output/graph/graph profiles ATCGATCGATCGATCG" << endl;
        cerr << endl;
        return 1;
    }
    
    string graph_path = argv[1];
    string profiles_dir = argv[2];
    string repeat_kmer = argv[3];
    
    cout << "Configuration:" << endl;
    cout << "  Graph path:    " << graph_path << endl;
    cout << "  Profiles dir:  " << profiles_dir << endl;
    cout << "  Repeat k-mer:  " << repeat_kmer << endl;
    cout << endl;
    
    // Load SDBG graph
    cout << "Loading SDBG graph..." << endl;
    SDBG sdbg;
    char* graph_cstr = new char[graph_path.length() + 1];
    strcpy(graph_cstr, graph_path.c_str());
    sdbg.LoadFromFile(graph_cstr);
    delete[] graph_cstr;
    
    cout << "✓ Graph loaded successfully!" << endl;
    cout << "  Graph size:    " << sdbg.size() << " nodes" << endl;
    cout << "  K-mer size:    " << sdbg.k() << endl;
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
    
    cout << "  Node multiplicity: " << sdbg.EdgeMultiplicity(repeat_node) << endl;
    cout << "  Node outdegree:    " << sdbg.EdgeOutdegree(repeat_node) << endl;
    cout << endl;
    
    // Create CAS workflow
    cout << "Initializing CAS workflow..." << endl;
    CasWorkflow workflow(sdbg, profiles_dir);
    
    cout << "✓ Workflow initialized with default parameters:" << endl;
    cout << "  FIRST_GENE_MIN_DIST:    50 bp" << endl;
    cout << "  FIRST_GENE_MAX_DIST:    1000 bp" << endl;
    cout << "  MAX_START_CANDIDATES:   5000" << endl;
    cout << "  MAX_LOCUS_BP:           41591 bp" << endl;
    cout << "  OVERLAP_ALLOWANCE:      15 bp" << endl;
    cout << "  INTERGENIC_MAX:         100 bp" << endl;
    cout << "  BEAM_WIDTH:             10" << endl;
    cout << "  MIN_NORMALIZED_SCORE:   0.0" << endl;
    cout << endl;
    
    // Run CAS gene detection
    cout << "============================================" << endl;
    cout << "Starting CAS gene detection..." << endl;
    cout << "============================================" << endl;
    cout << endl;
    
    auto detected_genes = workflow.DetectCasGenes(repeat_node);
    
    cout << endl;
    cout << "============================================" << endl;
    cout << "RESULTS" << endl;
    cout << "============================================" << endl;
    cout << endl;
    
    cout << "Total genes detected: " << detected_genes.size() << endl;
    
    if (detected_genes.empty()) {
        cout << endl;
        cout << "No CAS genes detected." << endl;
        cout << "This could mean:" << endl;
        cout << "  - No start codons found in search range" << endl;
        cout << "  - No genes scored above threshold" << endl;
        cout << "  - Wrong repeat node or direction" << endl;
        return 0;
    }
    
    // Calculate total operon length
    int total_operon_length = 0;
    if (!detected_genes.empty()) {
        const auto& last_gene = detected_genes.back();
        total_operon_length = last_gene.distance_from_repeat + last_gene.gene_length;
    }
    
    cout << "Total operon length: " << total_operon_length << " bp" << endl;
    cout << endl;
    
    cout << "Detected CAS Genes:" << endl;
    cout << "-------------------" << endl;
    cout << endl;
    
    for (size_t i = 0; i < detected_genes.size(); i++) {
        const auto& gene = detected_genes[i];
        
        cout << "Gene #" << (i + 1) << ":" << endl;
        cout << "  Profile:           " << gene.profile_name << endl;
        cout << "  Start node:        " << gene.start_node << endl;
        cout << "  End node:          " << gene.end_node << endl;
        cout << "  Distance from R:   " << gene.distance_from_repeat << " bp" << endl;
        cout << "  Gene length:       " << gene.gene_length << " bp" << endl;
        cout << "  Bit score:         " << fixed << setprecision(2) << gene.bit_score << endl;
        cout << "  Normalized score:  " << fixed << setprecision(4) << gene.normalized_score << endl;
        cout << "  Amino acids:       " << gene.amino_acids.length() << " AA" << endl;
        
        // Show amino acid sequence (first 30 AA)
        if (!gene.amino_acids.empty()) {
            cout << "  AA sequence:       ";
            int num_to_show = min(30, (int)gene.amino_acids.length());
            for (int j = 0; j < num_to_show; j++) {
                cout << gene.amino_acids[j];
            }
            if (gene.amino_acids.length() > 30) {
                cout << "... (" << gene.amino_acids.length() << " total)";
            }
            cout << endl;
        }
        
        cout << "  Node path length:  " << gene.node_path.size() << " nodes" << endl;
        cout << endl;
    }
    
    // Summary statistics
    cout << "============================================" << endl;
    cout << "SUMMARY" << endl;
    cout << "============================================" << endl;
    cout << endl;
    
    double avg_gene_length = 0.0;
    double avg_normalized_score = 0.0;
    int total_aa = 0;
    
    for (const auto& gene : detected_genes) {
        avg_gene_length += gene.gene_length;
        avg_normalized_score += gene.normalized_score;
        total_aa += gene.amino_acids.length();
    }
    
    if (!detected_genes.empty()) {
        avg_gene_length /= detected_genes.size();
        avg_normalized_score /= detected_genes.size();
    }
    
    cout << "Average gene length:      " << fixed << setprecision(1) << avg_gene_length << " bp" << endl;
    cout << "Average normalized score: " << fixed << setprecision(4) << avg_normalized_score << endl;
    cout << "Total amino acids:        " << total_aa << " AA" << endl;
    cout << "Average AA per gene:      " << fixed << setprecision(1) 
         << (detected_genes.empty() ? 0.0 : (double)total_aa / detected_genes.size()) << " AA" << endl;
    cout << endl;
    
    cout << "Test completed successfully!" << endl;
    
    return 0;
}
