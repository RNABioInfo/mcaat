#include "profile.h"
#include <iostream>
#include <vector>

int main() {
    Profile hmm;
    hmm.LoadFromFile("hmm_test.hmm");
    
    // Test sequence with stop codons
    std::vector<std::string> bad_seq = {"M", "K", "L", "*", "A", "G", "*", "T"};
    
    std::cout << "Testing sequence with stop codons (*): ";
    for (const auto& aa : bad_seq) std::cout << aa;
    std::cout << std::endl;
    
    auto [score, path, matches] = hmm.ViterbiAlign(bad_seq);
    
    std::cout << "Score: " << score << " bits" << std::endl;
    std::cout << "Path length: " << path.length() << std::endl;
    std::cout << "HMM matches: " << matches << std::endl;
    
    std::cout << "\nExpected: -inf or very large negative score (paths with stop codons filtered)" << std::endl;
    
    return 0;
}
