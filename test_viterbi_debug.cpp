#include "profile.h"
#include <iostream>
#include <vector>

int main() {
    // Load HMM
    Profile hmm;
    hmm.LoadFromFile("hmm_test.hmm");
    
    // Simple test sequence (first 20 AAs from the 407 AA sequence)
    std::vector<std::string> seq = {"M","K","G","I","V","V","G","S","T","G","V","G","Q","H","Y","T","V","R","K","K"};
    
    // Run Viterbi
    auto [score, path, matches] = hmm.ViterbiAlign(seq);
    
    std::cout << "Score: " << score << std::endl;
    std::cout << "Matches: " << matches << " / " << hmm.GetLength() << std::endl;
    std::cout << "Path: " << path << std::endl;
    
    return 0;
}
