#include "profile.h"
#include <iostream>
#include <vector>
#include <string>

int main() {
    // Load HMM
    Profile hmm;
    hmm.LoadFromFile("hmm_test.hmm");
    
    // Sequence from 15.fasta (COG1518-sample1)
    std::string seq_str = "AIQTQSNLLEDAITTVNVRGGNVHVKASMRRRCPVKQIDQIMLLGSPVIFTMVLMCVSKQELPLHFFENFGKFCGRLSPRVSMASAIALNEQCRAAFDAHGLRCSHNEVEGPVYHLQQANNKAEEYSIVFDQVRDSFGAVRVKFGNRLQVAAMAELEFAETSDKRNGEGQARTKCNQKISDQTDLDHPMFTEANTDSQEDTTNKTLSVLGSTDTGNLLDATVDLGYLFDEGFFHEGRELSFTLATDVAEIELFRSTAVDRTVRKHCNSLLTPNEAVGIEAAHLTEDHVTALSQPGVGGGVGGSADKLPMDFVSSEQVRDEERIKFERIRHKIPYNRLTNVEPGEIGHKEKLGAYDREPVKRT";
    
    // Convert to vector<string> (one AA per element)
    std::vector<std::string> seq;
    for (char c : seq_str) {
        seq.push_back(std::string(1, c));
    }
    
    std::cout << "Sequence length: " << seq.size() << std::endl;
    
    // Run Viterbi
    auto [score, path, matches] = hmm.ViterbiAlign(seq);
    
    std::cout << "Score: " << score << std::endl;
    std::cout << "HMM matches (Viterbi): " << matches << " / " << hmm.GetLength() << std::endl;
    
    // Count M states in path
    int m_count = 0;
    for (char c : path) {
        if (c == 'M') m_count++;
    }
    std::cout << "M states in path: " << m_count << std::endl;
    
    return 0;
}
