#include <iostream>
#include "profile.h"

int main(int argc, char** argv) {
    std::string hmm_file = "hmm_test.hmm";
    
    if (argc > 1) {
        hmm_file = argv[1];
    }
    
    std::cout << "Loading HMM profile from: " << hmm_file << std::endl;
    std::cout << "========================================" << std::endl;
    
    Profile profile;
    if (!profile.LoadFromFile(hmm_file)) {
        std::cerr << "Failed to load profile!" << std::endl;
        return 1;
    }
    
    std::cout << "\n";
    profile.PrintSummary();
    
    // Test some queries
    std::cout << "\n=== Testing Queries ===" << std::endl;
    
    if (profile.GetLength() > 0) {
        // Test emission scores
        std::cout << "\nMatch emission for position 1, amino acid 'L': " 
                  << profile.GetMatchEmission(1, 'L') << std::endl;
        std::cout << "Match emission for position 1, amino acid 'A': " 
                  << profile.GetMatchEmission(1, 'A') << std::endl;
        
        // Test transition scores
        std::cout << "\nTransition M->M at position 1: " 
                  << profile.GetTransition(1, 2, 'M', 'M') << std::endl;
        std::cout << "Transition M->I at position 1: " 
                  << profile.GetTransition(1, 1, 'M', 'I') << std::endl;
        std::cout << "Transition M->D at position 1: " 
                  << profile.GetTransition(1, 2, 'M', 'D') << std::endl;
        
        // Show consensus sequence
        std::cout << "\nConsensus sequence (first 50 positions): ";
        for (int i = 1; i <= std::min(50, profile.GetLength()); ++i) {
            std::cout << profile.GetState(i).consensus;
        }
        std::cout << std::endl;
    }
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Test completed successfully!" << std::endl;
    
    return 0;
}
