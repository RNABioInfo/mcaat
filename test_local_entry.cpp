#include "profile.h"
#include <iostream>
#include <vector>
#include <string>

int main() {
    Profile hmm;
    hmm.LoadFromFile("hmm_test.hmm");
    
    std::string seq_str = "AIQTQSNLLEDAITTVNVRGGNVHVKASMRRRCPVKQIDQIMLLGSPVIFTMVLMCVSKQELPLHFFENFGKFCGRLSPRVSMASAIALNEQCRAAFDAHGLRCSHNEVEGPVYHLQQANNKAEEYSIVFDQVRDSFGAVRVKFGNRLQVAAMAELEFAETSDKRNGEGQARTKCNQKISDQTDLDHPMFTEANTDSQEDTTNKTLSVLGSTDTGNLLDATVDLGYLFDEGFFHEGRELSFTLATDVAEIELFRSTAVDRTVRKHCNSLLTPNEAVGIEAAHLTEDHVTALSQPGVGGGVGGSADKLPMDFVSSEQVRDEERIKFERIRHKIPYNRLTNVEPGEIGHKEKLGAYDREPVKRT";
    
    std::vector<std::string> seq;
    for (char c : seq_str) {
        seq.push_back(std::string(1, c));
    }
    
    auto [score, path, matches] = hmm.ViterbiAlign(seq);
    
    // Find where we entered and exited the HMM
    int entry_j = -1, exit_j = -1;
    for (size_t i = 0; i < path.size(); i++) {
        if (path[i] == 'M' && entry_j == -1) {
            // First M state - count how many states before this
            int j = 1;
            for (size_t k = 0; k < i; k++) {
                if (path[k] == 'M' || path[k] == 'D') j++;
            }
            entry_j = j;
        }
        if (path[path.size() - 1 - i] == 'M' && exit_j == -1) {
            // Last M state from the end
            int j = matches;
            for (size_t k = path.size() - 1 - i + 1; k < path.size(); k++) {
                if (path[k] == 'M' || path[k] == 'D') j++;
            }
            exit_j = j;
        }
    }
    
    std::cout << "Score: " << score << " bits" << std::endl;
    std::cout << "HMM coverage: " << entry_j << " - " << exit_j << " (" << (exit_j - entry_j + 1) << " positions)" << std::endl;
    std::cout << "M states: " << matches << std::endl;
    std::cout << "HMMER was: 13 - 308 (296 positions), score 105.5 bits" << std::endl;
    
    return 0;
}
