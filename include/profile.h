#ifndef PROFILE_H
#define PROFILE_H

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

// HMMER3 profile HMM structure
struct ProfileState {
    int position;  // Match state position
    std::vector<double> match_emissions;  // 20 amino acids
    std::vector<double> insert_emissions; // 20 amino acids (second line)
    std::vector<double> transitions;      // 7 transition probabilities (third line)
    char consensus;  // Consensus amino acid
};

class Profile {
private:
    std::string name_;
    int length_;  // Number of match states
    std::string alphabet_;  // "amino" or "DNA"
    
    // Header metadata
    std::map<std::string, std::string> metadata_;
    
    // Amino acid alphabet order (for indexing)
    std::vector<char> aa_alphabet_;
    
    // Profile states
    std::vector<ProfileState> states_;
    
    // COMPO line (composition/background frequencies)
    std::vector<double> compo_match_;
    std::vector<double> compo_insert_;
    std::vector<double> compo_transitions_;
    
    // Parse helper functions
    bool ParseHeader(const std::string& line);
    bool ParseHMMLine(const std::string& line);
    bool ParseStateLine(const std::vector<std::string>& lines, int& line_idx);
    std::vector<double> ParseScoreLine(const std::string& line);
    std::vector<std::string> SplitString(const std::string& str);

public:
    Profile();
    
    // Load from HMMER3 file
    bool LoadFromFile(const std::string& filename);
    
    // Getters
    std::string GetName() const { return name_; }
    int GetLength() const { return length_; }
    std::string GetAlphabet() const { return alphabet_; }
    
    const std::vector<ProfileState>& GetStates() const { return states_; }
    const ProfileState& GetState(int position) const;
    
    // Get emission score for amino acid at position
    double GetMatchEmission(int position, char amino_acid) const;
    double GetInsertEmission(int position, char amino_acid) const;
    
    // Get log-odds emission score (relative to null model)
    double GetMatchLogOdds(int position, char amino_acid) const;
    
    // Get bit score from log-odds
    double LogOddsToBits(double log_odds) const;
    
    // Viterbi alignment of amino acid sequence to HMM
    // Returns: <bit_score, best_path_states, hmm_end_position>
    std::tuple<double, std::string, int> ViterbiAlign(const std::vector<std::string>& aa_sequence) const;
    
    // Get transition score
    double GetTransition(int from_state, int to_state, char from_type, char to_type) const;
    
    // Print summary
    void PrintSummary() const;
    
    // Amino acid to index mapping
    int AminoAcidToIndex(char aa) const;
};

#endif // PROFILE_H
