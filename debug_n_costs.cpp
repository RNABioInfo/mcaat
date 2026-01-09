#include "profile.h"
#include <iostream>

int main() {
    Profile hmm;
    hmm.LoadFromFile("hmm_test.hmm");
    
    // Check what the N state costs are
    std::cout << "N->N transition: " << hmm.GetSpecialTransition(0, 1) << std::endl;
    std::cout << "N->B transition: " << hmm.GetSpecialTransition(0, 0) << std::endl;
    
    // Simulate N state accumulation for 13 positions
    std::string seq = "AIQTQSNLLEDAITTVNVRGGNVHVKASMRRRCPVKQIDQIMLLGSPVIFTMVLMCVSKQELPLHFFENFGKFCGRLSPRVSMASAIALNEQCRAAFDAHGLRCSHNEVEGPVYHLQQANNKAEEYSIVFDQVRDSFGAVRVKFGNRLQVAAMAELEFAETSDKRNGEGQARTKCNQKISDQTDLDHPMFTEANTDSQEDTTNKTLSVLGSTDTGNLLDATVDLGYLFDEGFFHEGRELSFTLATDVAEIELFRSTAVDRTVRKHCNSLLTPNEAVGIEAAHLTEDHVTALSQPGVGGGVGGSADKLPMDFVSSEQVRDEERIKFERIRHKIPYNRLTNVEPGEIGHKEKLGAYDREPVKRT";
    
    return 0;
}
