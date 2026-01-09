#include <iostream>
#include <cmath>

int main() {
    std::cout << "=== Exit Penalty Analysis ===" << std::endl;
    std::cout << "Formula: log(2.0 / (length + 2)) * 2" << std::endl;
    std::cout << std::endl;
    
    for (int len : {6, 10, 50, 100, 362, 1000}) {
        double penalty_nats = std::log(2.0 / (len + 2)) * 2.0;
        double penalty_bits = penalty_nats / std::log(2.0);
        
        std::cout << "Length " << len << ":" << std::endl;
        std::cout << "  Penalty: " << penalty_bits << " bits" << std::endl;
        std::cout << "  Effect: ";
        if (penalty_bits < 0) {
            std::cout << "DECREASES score (penalty)" << std::endl;
        } else {
            std::cout << "INCREASES score (bonus)" << std::endl;
        }
    }
    
    std::cout << std::endl;
    std::cout << "Observation: Longer sequences get MORE NEGATIVE penalty" << std::endl;
    std::cout << "This means: We PENALIZE long sequences more than short ones" << std::endl;
    std::cout << std::endl;
    std::cout << "Problem: Raw HMM scores are ALREADY more negative for longer sequences!" << std::endl;
    std::cout << "So we're DOUBLE-PENALIZING length." << std::endl;
    
    return 0;
}
