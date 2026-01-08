#include <iostream>
#include <string>
#include <algorithm>
#include <cstring>
#include <vector>
#include "sdbg/sdbg.h"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_prefix> <kmer_sequence>" << std::endl;
        return 1;
    }

    std::string graph_prefix = argv[1];
    std::string kmer = argv[2];

    // Load the graph
    std::cout << "Loading graph from: " << graph_prefix << std::endl;
    SDBG sdbg;
    char* cstr = new char[graph_prefix.length() + 1];
    strcpy(cstr, graph_prefix.c_str());
    sdbg.LoadFromFile(cstr);
    delete[] cstr;
    
    std::cout << "Graph loaded. K-mer size: " << sdbg.k() << std::endl;
    std::cout << "Number of nodes: " << sdbg.size() << std::endl;
    
    if (kmer.length() != (size_t)sdbg.k()) {
        std::cerr << "Error: K-mer length (" << kmer.length() 
                  << ") doesn't match graph k (" << sdbg.k() << ")" << std::endl;
        return 1;
    }

    // Search for the k-mer
    std::cout << "\nSearching for k-mer: " << kmer << std::endl;
    
    // Convert DNA sequence to encoded format (A=1, C=2, G=3, T=4)
    uint8_t *encoded = new uint8_t[sdbg.k()];
    for (uint32_t i = 0; i < sdbg.k(); i++) {
        switch(kmer[i]) {
            case 'A': case 'a': encoded[i] = 1; break;
            case 'C': case 'c': encoded[i] = 2; break;
            case 'G': case 'g': encoded[i] = 3; break;
            case 'T': case 't': encoded[i] = 4; break;
            default: 
                std::cerr << "Invalid base: " << kmer[i] << std::endl;
                delete[] encoded;
                return 1;
        }
    }
    
    // Try to find the k-mer
    uint64_t node_id = sdbg.IndexBinarySearch(encoded);
    
    if (node_id != SDBG::kNullID) {
        std::cout << "✓ Found! Node ID: " << node_id << std::endl;
        
        // Verify by reading back the k-mer
        std::vector<uint8_t> seq(sdbg.k());
        sdbg.GetLabel(node_id, seq.data());
        std::string read_seq;
        for (uint32_t i = 0; i < sdbg.k(); i++) {
            read_seq += "ACGT"[seq[i] - 1];
        }
        std::cout << "  Verification - Node " << node_id << " k-mer: " << read_seq << std::endl;
    } else {
        std::cout << "✗ K-mer not found in graph." << std::endl;
        
        // Try reverse complement
        std::string rc_kmer = kmer;
        for (size_t i = 0; i < rc_kmer.length(); i++) {
            switch(rc_kmer[i]) {
                case 'A': case 'a': rc_kmer[i] = 'T'; break;
                case 'T': case 't': rc_kmer[i] = 'A'; break;
                case 'C': case 'c': rc_kmer[i] = 'G'; break;
                case 'G': case 'g': rc_kmer[i] = 'C'; break;
            }
        }
        std::reverse(rc_kmer.begin(), rc_kmer.end());
        
        std::cout << "\nTrying reverse complement: " << rc_kmer << std::endl;
        
        // Encode RC
        for (uint32_t i = 0; i < sdbg.k(); i++) {
            switch(rc_kmer[i]) {
                case 'A': case 'a': encoded[i] = 1; break;
                case 'C': case 'c': encoded[i] = 2; break;
                case 'G': case 'g': encoded[i] = 3; break;
                case 'T': case 't': encoded[i] = 4; break;
            }
        }
        
        node_id = sdbg.IndexBinarySearch(encoded);
        
        if (node_id != SDBG::kNullID) {
            std::cout << "✓ Found RC! Node ID: " << node_id << std::endl;
            std::vector<uint8_t> seq(sdbg.k());
            sdbg.GetLabel(node_id, seq.data());
            std::string read_seq;
            for (uint32_t i = 0; i < sdbg.k(); i++) {
                read_seq += "ACGT"[seq[i] - 1];
            }
            std::cout << "  Verification - Node " << node_id << " k-mer: " << read_seq << std::endl;
        } else {
            std::cout << "✗ Reverse complement also not found." << std::endl;
        }
    }
    
    delete[] encoded;
    
    return 0;
}
