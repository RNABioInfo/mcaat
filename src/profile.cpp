#include "profile.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cmath>

Profile::Profile() : length_(0), alphabet_("amino") {
    // Standard HMMER amino acid order
    aa_alphabet_ = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
}

std::vector<std::string> Profile::SplitString(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (iss >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<double> Profile::ParseScoreLine(const std::string& line) {
    std::vector<double> scores;
    std::istringstream iss(line);
    std::string token;
    
    while (iss >> token) {
        if (token == "*") {
            scores.push_back(INFINITY);  // Special end marker
        } else {
            try {
                // HMMER3 format stores positive values representing -log(P)
                // Convert to log(P) by negating
                double value = std::stod(token);
                scores.push_back(-value);  // Negate to get log probability
            } catch (...) {
                // Skip non-numeric tokens
            }
        }
    }
    return scores;
}

bool Profile::ParseHeader(const std::string& line) {
    auto tokens = SplitString(line);
    if (tokens.empty()) return false;
    
    if (tokens[0] == "NAME" && tokens.size() > 1) {
        name_ = tokens[1];
    } else if (tokens[0] == "LENG" && tokens.size() > 1) {
        length_ = std::stoi(tokens[1]);
    } else if (tokens[0] == "ALPH" && tokens.size() > 1) {
        alphabet_ = tokens[1];
    } else {
        // Store other metadata
        if (tokens.size() > 1) {
            metadata_[tokens[0]] = tokens[1];
        }
    }
    return true;
}

bool Profile::LoadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }
    
    std::string line;
    std::vector<std::string> file_lines;
    
    // Read all lines
    while (std::getline(file, line)) {
        file_lines.push_back(line);
    }
    file.close();
    
    if (file_lines.empty()) {
        std::cerr << "Empty file" << std::endl;
        return false;
    }
    
    // Parse header
    int line_idx = 0;
    while (line_idx < file_lines.size()) {
        const std::string& current_line = file_lines[line_idx];
        
        if (current_line.find("HMM") == 0 && current_line.find("A") != std::string::npos) {
            // Found HMM header line with amino acids
            line_idx++;
            break;
        }
        
        ParseHeader(current_line);
        line_idx++;
    }
    
    // Skip transition header line (m->m m->i m->d ...)
    if (line_idx < file_lines.size()) {
        line_idx++;
    }
    
    // Parse COMPO line
    if (line_idx < file_lines.size() && file_lines[line_idx].find("COMPO") != std::string::npos) {
        auto tokens = SplitString(file_lines[line_idx]);
        if (tokens.size() > 1) {
            // HMMER3 format stores positive values, negate to get log-probabilities
            for (size_t i = 1; i < tokens.size() && i <= 20; ++i) {
                double value = std::stod(tokens[i]);
                compo_match_.push_back(-value);  // Negate to get log-probability
            }
        }
        line_idx++;
        
        // Next line is insert emissions for COMPO
        if (line_idx < file_lines.size()) {
            auto scores = ParseScoreLine(file_lines[line_idx]);
            if (scores.size() >= 20) {
                compo_insert_.assign(scores.begin(), scores.begin() + 20);
            }
            line_idx++;
        }
        
        // Next line is transitions for COMPO
        if (line_idx < file_lines.size()) {
            compo_transitions_ = ParseScoreLine(file_lines[line_idx]);
            line_idx++;
        }
    }
    
    // Parse states
    while (line_idx < file_lines.size()) {
        const std::string& current_line = file_lines[line_idx];
        
        // Check for end marker
        if (current_line.find("//") == 0) {
            break;
        }
        
        // Skip empty lines
        if (current_line.empty() || current_line.find_first_not_of(" \t") == std::string::npos) {
            line_idx++;
            continue;
        }
        
        // Try to parse a state (3 lines per state)
        auto tokens = SplitString(current_line);
        
        // First line should start with a position number
        if (!tokens.empty() && std::isdigit(tokens[0][0])) {
            ProfileState state;
            
            try {
                state.position = std::stoi(tokens[0]);
                
                // Parse match emissions (20 amino acids)
                // HMMER3 format stores positive values representing -log(P)
                // Negate to get log(P)
                for (size_t i = 1; i < tokens.size() && i <= 20; ++i) {
                    double value = std::stod(tokens[i]);
                    state.match_emissions.push_back(-value);  // Negate to get log-probability
                }
                
                // Get consensus amino acid - it's after position number (22nd token)
                // Format: "position AA1 AA2 ... AA20 position consensus - - -"
                if (tokens.size() > 22) {
                    state.consensus = tokens[22][0];
                } else if (tokens.size() > 21) {
                    state.consensus = tokens[21][0];
                }
                
                line_idx++;
                
                // Parse insert emissions (second line)
                if (line_idx < file_lines.size()) {
                    auto insert_scores = ParseScoreLine(file_lines[line_idx]);
                    if (insert_scores.size() >= 20) {
                        state.insert_emissions.assign(insert_scores.begin(), insert_scores.begin() + 20);
                    }
                    line_idx++;
                }
                
                // Parse transitions (third line)
                if (line_idx < file_lines.size()) {
                    state.transitions = ParseScoreLine(file_lines[line_idx]);
                    line_idx++;
                }
                
                states_.push_back(state);
                
            } catch (const std::exception& e) {
                std::cerr << "Error parsing state at line " << line_idx << ": " << e.what() << std::endl;
                line_idx++;
            }
        } else {
            line_idx++;
        }
    }
    
    std::cout << "Loaded profile: " << name_ << std::endl;
    std::cout << "Length: " << length_ << " states" << std::endl;
    std::cout << "Parsed: " << states_.size() << " states" << std::endl;
    std::cout << "Alphabet: " << alphabet_ << std::endl;
    
    return !states_.empty();
}

int Profile::AminoAcidToIndex(char aa) const {
    aa = std::toupper(aa);
    auto it = std::find(aa_alphabet_.begin(), aa_alphabet_.end(), aa);
    if (it != aa_alphabet_.end()) {
        return std::distance(aa_alphabet_.begin(), it);
    }
    return -1;  // Invalid amino acid
}

const ProfileState& Profile::GetState(int position) const {
    if (position < 1 || position > static_cast<int>(states_.size())) {
        throw std::out_of_range("Invalid state position");
    }
    return states_[position - 1];  // States are 1-indexed
}

double Profile::GetMatchEmission(int position, char amino_acid) const {
    const auto& state = GetState(position);
    int idx = AminoAcidToIndex(amino_acid);
    if (idx < 0 || idx >= static_cast<int>(state.match_emissions.size())) {
        return INFINITY;  // Invalid amino acid
    }
    return state.match_emissions[idx];
}

double Profile::GetInsertEmission(int position, char amino_acid) const {
    const auto& state = GetState(position);
    int idx = AminoAcidToIndex(amino_acid);
    if (idx < 0 || idx >= static_cast<int>(state.insert_emissions.size())) {
        return -INFINITY;  // Invalid amino acid
    }
    // Return raw log-probability (already negated from HMMER3 format)
    return state.insert_emissions[idx];
}

double Profile::GetMatchLogOdds(int position, char amino_acid) const {
    const auto& state = GetState(position);
    int idx = AminoAcidToIndex(amino_acid);
    if (idx < 0 || idx >= static_cast<int>(state.match_emissions.size())) {
        return -INFINITY;  // Invalid amino acid
    }
    
    // Get emission score (log probability)
    double emission = state.match_emissions[idx];
    
    // Get background frequency (null model)
    double background = 0.0;
    if (idx < static_cast<int>(compo_match_.size())) {
        background = compo_match_[idx];
    } else {
        // Uniform background if not available
        background = std::log(1.0 / 20.0);
    }
    
    // Log-odds = log(P_match) - log(P_background) = log(P_match / P_background)
    return emission - background;
}

double Profile::LogOddsToBits(double log_odds) const {
    // HMMER bit scores: bits = log-odds / log(2)
    return log_odds / std::log(2.0);
}

double Profile::GetTransition(int from_state, int to_state, char from_type, char to_type) const {
    // Transition indices in HMMER3:
    // 0: m->m, 1: m->i, 2: m->d, 3: i->m, 4: i->i, 5: d->m, 6: d->d
    const auto& state = GetState(from_state);
    
    int trans_idx = -1;
    if (from_type == 'm' || from_type == 'M') {
        if (to_type == 'm' || to_type == 'M') trans_idx = 0;
        else if (to_type == 'i' || to_type == 'I') trans_idx = 1;
        else if (to_type == 'd' || to_type == 'D') trans_idx = 2;
    } else if (from_type == 'i' || from_type == 'I') {
        if (to_type == 'm' || to_type == 'M') trans_idx = 3;
        else if (to_type == 'i' || to_type == 'I') trans_idx = 4;
    } else if (from_type == 'd' || from_type == 'D') {
        if (to_type == 'm' || to_type == 'M') trans_idx = 5;
        else if (to_type == 'd' || to_type == 'D') trans_idx = 6;
    }
    
    if (trans_idx < 0 || trans_idx >= static_cast<int>(state.transitions.size())) {
        return INFINITY;
    }
    
    return state.transitions[trans_idx];
}

std::tuple<double, std::string, int> Profile::ViterbiAlign(const std::vector<std::string>& aa_sequence) const {
    int seq_len = aa_sequence.size();
    int hmm_len = length_;
    
    // DP matrix: [sequence_position][hmm_position][state_type]
    // state_type: 0=Match, 1=Insert, 2=Delete
    const int NUM_STATES = 3;
    const int M = 0, I = 1, D = 2;
    
    // Initialize DP table with -infinity (log space)
    std::vector<std::vector<std::vector<double>>> dp(
        seq_len + 1, 
        std::vector<std::vector<double>>(hmm_len + 1, std::vector<double>(NUM_STATES, -INFINITY))
    );
    
    // Traceback matrix
    std::vector<std::vector<std::vector<int>>> traceback(
        seq_len + 1,
        std::vector<std::vector<int>>(hmm_len + 1, std::vector<int>(NUM_STATES, -1))
    );
    
    // Initialize: start at position 0 in match state
    dp[0][0][M] = 0.0;
    
    // Fill DP matrix
    for (int i = 0; i <= seq_len; ++i) {
        for (int j = 0; j <= hmm_len; ++j) {
            // Match state: consume AA and advance HMM
            if (i > 0 && j > 0 && j <= hmm_len) {
                char aa = aa_sequence[i-1][0];  // Get amino acid
                double emission = GetMatchEmission(j, aa);  // Use raw log-probability
                
                // From previous Match (j-1 must be >= 1 for GetTransition)
                double from_m = -INFINITY;
                if (j > 1) {
                    from_m = dp[i-1][j-1][M] + emission + GetTransition(j-1, j, 'M', 'M');
                } else {
                    // First position (j=1), come from initial state
                    from_m = dp[i-1][j-1][M] + emission;
                }
                
                // From previous Insert
                double from_i = -INFINITY;
                if (j > 1) {
                    from_i = dp[i-1][j-1][I] + emission + GetTransition(j-1, j, 'I', 'M');
                }
                
                // From previous Delete
                double from_d = -INFINITY;
                if (j > 1) {
                    from_d = dp[i-1][j-1][D] + emission + GetTransition(j-1, j, 'D', 'M');
                }
                
                double best = std::max({from_m, from_i, from_d});
                if (best > dp[i][j][M]) {
                    dp[i][j][M] = best;
                    if (best == from_m) traceback[i][j][M] = M;
                    else if (best == from_i) traceback[i][j][M] = I;
                    else traceback[i][j][M] = D;
                }
            }
            
            // Insert state: consume AA but don't advance HMM
            if (i > 0 && j > 0 && j <= hmm_len) {
                char aa = aa_sequence[i-1][0];
                // Get raw insert emission (log-probability)
                const auto& state = GetState(j);
                int idx = AminoAcidToIndex(aa);
                double emission = (idx >= 0 && idx < static_cast<int>(state.insert_emissions.size())) 
                    ? state.insert_emissions[idx] : -INFINITY;
                
                // From previous Match
                double from_m = dp[i-1][j][M] + emission + GetTransition(j, j, 'M', 'I');
                // From previous Insert
                double from_i = dp[i-1][j][I] + emission + GetTransition(j, j, 'I', 'I');
                
                double best = std::max(from_m, from_i);
                if (best > dp[i][j][I]) {
                    dp[i][j][I] = best;
                    traceback[i][j][I] = (best == from_m) ? M : I;
                }
            }
            
            // Delete state: advance HMM but don't consume AA
            if (j > 1 && i <= seq_len && j <= hmm_len) {
                // From previous Match (j-1 >= 1)
                double from_m = dp[i][j-1][M] + GetTransition(j-1, j, 'M', 'D');
                // From previous Delete
                double from_d = dp[i][j-1][D] + GetTransition(j-1, j, 'D', 'D');
                
                double best = std::max(from_m, from_d);
                if (best > dp[i][j][D]) {
                    dp[i][j][D] = best;
                    traceback[i][j][D] = (best == from_m) ? M : D;
                }
            }
        }
    }
    
    // Find best final score (local alignment - can end anywhere)
    double best_score = -INFINITY;
    int best_final_state = M;
    int best_final_pos = hmm_len;
    
    // Check all possible ending positions and states
    for (int j = 1; j <= hmm_len; ++j) {
        if (dp[seq_len][j][M] > best_score) {
            best_score = dp[seq_len][j][M];
            best_final_state = M;
            best_final_pos = j;
        }
        if (dp[seq_len][j][I] > best_score) {
            best_score = dp[seq_len][j][I];
            best_final_state = I;
            best_final_pos = j;
        }
        if (dp[seq_len][j][D] > best_score) {
            best_score = dp[seq_len][j][D];
            best_final_state = D;
            best_final_pos = j;
        }
    }
    
    // Traceback to get path
    std::string path;
    int i = seq_len, j = best_final_pos, state = best_final_state;
    while (i > 0 || j > 0) {
        if (state == M) {
            path = "M" + path;
            int prev = traceback[i][j][M];
            i--; j--;
            state = prev;
        } else if (state == I) {
            path = "I" + path;
            int prev = traceback[i][j][I];
            i--;
            state = prev;
        } else if (state == D) {
            path = "D" + path;
            int prev = traceback[i][j][D];
            j--;
            state = prev;
        } else {
            break;
        }
    }
    
    // Convert log-probability score to bits with null model correction
    // Null model score: what would a random sequence score?
    double null_score = 0.0;
    for (const auto& aa : aa_sequence) {
        int idx = AminoAcidToIndex(aa[0]);
        if (idx >= 0 && idx < static_cast<int>(compo_match_.size())) {
            null_score += compo_match_[idx];  // Background frequency (log-prob)
        }
    }
    
    // Log-odds = HMM score - null model score
    double log_odds = best_score - null_score;
    
    // Convert to bits (HMMER-compatible)
    double bit_score = log_odds / std::log(2.0);
    
    return {bit_score, path, best_final_pos};
}

void Profile::PrintSummary() const {
    std::cout << "=== Profile HMM Summary ===" << std::endl;
    std::cout << "Name: " << name_ << std::endl;
    std::cout << "Length: " << length_ << std::endl;
    std::cout << "Alphabet: " << alphabet_ << std::endl;
    std::cout << "States loaded: " << states_.size() << std::endl;
    
    if (!states_.empty()) {
        std::cout << "\nFirst state:" << std::endl;
        const auto& first = states_[0];
        std::cout << "  Position: " << first.position << std::endl;
        std::cout << "  Consensus: " << first.consensus << std::endl;
        std::cout << "  Match emissions (first 5): ";
        for (size_t i = 0; i < std::min(size_t(5), first.match_emissions.size()); ++i) {
            std::cout << first.match_emissions[i] << " ";
        }
        std::cout << std::endl;
        
        std::cout << "\nLast state:" << std::endl;
        const auto& last = states_.back();
        std::cout << "  Position: " << last.position << std::endl;
        std::cout << "  Consensus: " << last.consensus << std::endl;
        std::cout << "  Match emissions (first 5): ";
        for (size_t i = 0; i < std::min(size_t(5), last.match_emissions.size()); ++i) {
            std::cout << last.match_emissions[i] << " ";
        }
        std::cout << std::endl;
    }
}
