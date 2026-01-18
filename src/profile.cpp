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
            scores.push_back(-INFINITY);  // Impossible transition in log-space
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
    
    // Initialize special state transitions (HMMER3 style)
    // Special states: N=0, B=1, E=2, C=3, J=4
    // Each has [MOVE, LOOP] transitions
    special_transitions_.resize(5, std::vector<double>(2, -INFINITY));
    
    // Default values for SINGLE-DOMAIN alignment (no multi-domain)
    // These are in log-probability space (negative values)
    special_transitions_[0][0] = -0.5;  // N->B (MOVE)
    special_transitions_[0][1] = -2.0;  // N->N (LOOP)
    special_transitions_[1][0] = 0.0;   // B->M (implicit, handled separately)
    special_transitions_[1][1] = -INFINITY;
    special_transitions_[2][0] = -0.5;  // E->C (MOVE)
    special_transitions_[2][1] = -INFINITY;  // E->J (LOOP) - DISABLED for single-domain!
    special_transitions_[3][0] = 0.0;   // C->T (MOVE) - FREE transition to terminal
    special_transitions_[3][1] = -2.0;  // C->C (LOOP)
    special_transitions_[4][0] = -INFINITY;  // J->B (MOVE) - DISABLED for single-domain!
    special_transitions_[4][1] = -INFINITY;  // J->J (LOOP) - DISABLED for single-domain!
    
    // Detect local vs glocal mode from STATS line
    is_local_ = true;  // Use local mode for flexible entry/exit
    
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
        // Invalid amino acid (e.g., stop codon *)
        // MegaGTA approach: filter out stop codons completely
        return -INFINITY;
    }
    return state.match_emissions[idx];
}

double Profile::GetInsertEmission(int position, char amino_acid) const {
    const auto& state = GetState(position);
    int idx = AminoAcidToIndex(amino_acid);
    if (idx < 0 || idx >= static_cast<int>(state.insert_emissions.size())) {
        // Invalid amino acid (e.g., stop codon *)
        // MegaGTA approach: filter out stop codons completely
        return -INFINITY;
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

double Profile::GetInsertLogOdds(int position, char amino_acid) const {
    const auto& state = GetState(position);
    int idx = AminoAcidToIndex(amino_acid);
    if (idx < 0 || idx >= static_cast<int>(state.insert_emissions.size())) {
        return -INFINITY;
    }
    
    double emission = state.insert_emissions[idx];
    
    // Background from insert composition or uniform
    double background = 0.0;
    if (idx < static_cast<int>(compo_insert_.size())) {
        background = compo_insert_[idx];
    } else {
        background = std::log(1.0 / 20.0);
    }
    
    return emission - background;
}

double Profile::LogOddsToBits(double log_odds) const {
    // HMMER bit scores: bits = log-odds / log(2)
    return log_odds / std::log(2.0);
}

double Profile::GetSpecialTransition(int state, int type) const {
    if (state < 0 || state >= static_cast<int>(special_transitions_.size()) ||
        type < 0 || type >= 2) {
        return -INFINITY;
    }
    return special_transitions_[state][type];
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
        return -INFINITY;
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
    
    // Special states: N=0, B=1, E=2, C=3, J=4
    const int N = 0, B = 1, E = 2, C = 3, J = 4;
    
    // For traceback, use negative values to mark special state entries
    const int FROM_B = -1;  // Entered from B state
    
    // Initialize DP table with -infinity (log space)
    std::vector<std::vector<std::vector<double>>> dp(
        seq_len + 1, 
        std::vector<std::vector<double>>(hmm_len + 1, std::vector<double>(NUM_STATES, -INFINITY))
    );
    
    // Special state matrix: [sequence_position][special_state]
    std::vector<std::vector<double>> xmx(
        seq_len + 1,
        std::vector<double>(5, -INFINITY)
    );
    
    // Traceback matrix
    std::vector<std::vector<std::vector<int>>> traceback(
        seq_len + 1,
        std::vector<std::vector<int>>(hmm_len + 1, std::vector<int>(NUM_STATES, -1))
    );
    
    // HMMER-style initialization (row 0)
    xmx[0][N] = 0.0;                              // S->N, p=1
    xmx[0][B] = xmx[0][N] + GetSpecialTransition(N, 0);  // S->N->B
    xmx[0][E] = xmx[0][C] = xmx[0][J] = -INFINITY;  // Need sequence
    
    // Exit score: 0 for local mode, -inf for glocal
    double esc = is_local_ ? 0.0 : -INFINITY;
    
    // Fill DP matrix
    for (int i = 1; i <= seq_len; ++i) {
        // Initialize row
        dp[i][0][M] = dp[i][0][I] = dp[i][0][D] = -INFINITY;
        xmx[i][E] = -INFINITY;
        
        for (int j = 1; j <= hmm_len; ++j) {
            // Match state: consume AA and advance HMM
            char aa = aa_sequence[i-1][0];  // Get amino acid
            double emission = GetMatchEmission(j, aa);  // Use raw log-probability
            
            // From previous Match
            double from_m = -INFINITY;
            if (j > 1) {
                from_m = dp[i-1][j-1][M] + emission + GetTransition(j-1, j, 'M', 'M');
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
            
            // From B state (entry): HMMER local alignment  
            // In local mode, B can enter at ANY match state
            // Entry choice is driven by match quality and N-state costs
            double from_b = xmx[i-1][B] + emission;  // B→Mj entry
            
            double best = std::max({from_m, from_i, from_d, from_b});
            dp[i][j][M] = best;
            if (best == from_m) traceback[i][j][M] = M;
            else if (best == from_i) traceback[i][j][M] = I;
            else if (best == from_d) traceback[i][j][M] = D;
            else traceback[i][j][M] = -2;  // FROM_B marker
            
            // E state update: can exit from any match state (local mode)
            xmx[i][E] = std::max(xmx[i][E], dp[i][j][M] + esc);
            
            // Insert state: consume AA but don't advance HMM
            if (i > 0 && j > 0 && j <= hmm_len) {
                char aa = aa_sequence[i-1][0];
                // Get raw insert emission (log-probability)
                const auto& state = GetState(j);
                int idx = AminoAcidToIndex(aa);
                double emission = (idx >= 0 && idx < static_cast<int>(state.insert_emissions.size())) 
                    ? state.insert_emissions[idx] : -100.0;  // Large penalty for invalid AA
                
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
            if (j > 1) {
                // From previous Match
                double from_m = dp[i][j-1][M] + GetTransition(j-1, j, 'M', 'D');
                // From previous Delete
                double from_d = dp[i][j-1][D] + GetTransition(j-1, j, 'D', 'D');
                
                double best = std::max(from_m, from_d);
                dp[i][j][D] = best;
                traceback[i][j][D] = (best == from_m) ? M : D;
            }
        }
        
        // Special state transitions (HMMER-style)
        // These happen after completing all HMM positions for this sequence position
        
        // J state: E->J or J->J
        double sc_j = std::max(
            xmx[i-1][J] + GetSpecialTransition(J, 1),  // J->J (LOOP)
            xmx[i][E] + GetSpecialTransition(E, 1)     // E->J (LOOP)
        );
        xmx[i][J] = sc_j;
        
        // N state: N->N with emission (N-terminal unaligned region)
        char aa = aa_sequence[i-1][0];
        int aa_idx = AminoAcidToIndex(aa);
        double null_emission = (aa_idx >= 0 && aa_idx < static_cast<int>(compo_match_.size())) 
            ? compo_match_[aa_idx] : -INFINITY;  // Filter stop codons (MegaGTA approach)
        
        xmx[i][N] = xmx[i-1][N] + GetSpecialTransition(N, 1) + null_emission;  // N->N (LOOP)
        
        // C state: E->C or C->C with emission (C-terminal unaligned region)
        double sc_c = std::max(
            xmx[i-1][C] + GetSpecialTransition(C, 1) + null_emission,  // C->C (LOOP)
            xmx[i][E] + GetSpecialTransition(E, 0)     // E->C (MOVE, no emission)
        );
        xmx[i][C] = sc_c;
        
        // B state: N->B or J->B
        double sc_b = std::max(
            xmx[i][N] + GetSpecialTransition(N, 0),  // N->B (MOVE)
            xmx[i][J] + GetSpecialTransition(J, 0)   // J->B (MOVE)
        );
        xmx[i][B] = sc_b;
    }
    
    // HMMER local alignment for beam search scoring:
    // We consume the full sequence, but allow flexible HMM alignment within it
    // Path: N*(0+) → B → M+ → E → C*(0+) where N+M+C consume all seq_len AAs
    
    double best_score = xmx[seq_len][C] + GetSpecialTransition(C, 0);  // Final score
    
    // For traceback: find best M→E at end of sequence
    // (The C state path from there to seq_len is deterministic)
    int best_final_i = seq_len;
    int best_final_pos = 1;
    int best_final_state = M;
    double best_e_score = -INFINITY;
    
    for (int j = 1; j <= hmm_len; ++j) {
        double e_score = dp[seq_len][j][M] + esc;
        if (e_score > best_e_score) {
            best_e_score = e_score;
            best_final_pos = j;
        }
    }
    
    // Traceback to get path
    std::string path;
    int i = best_final_i, j = best_final_pos;
    int state = best_final_state;
    
    while (i > 0 && j > 0) {
        if (state == M) {
            path = "M" + path;
            int prev = traceback[i][j][M];
            if (prev == -2) {
                // FROM_B: Entered HMM from B state
                break;
            }
            i--; j--;
            if (prev < 0) {
                break;  // Invalid/uninitialized
            }
            state = prev;
        } else if (state == I) {
            path = "I" + path;
            int prev = traceback[i][j][I];
            i--;
            if (prev < 0) break;
            state = prev;
        } else if (state == D) {
            path = "D" + path;
            int prev = traceback[i][j][D];
            j--;
            if (prev < 0) break;
            state = prev;
        } else {
            break;
        }
    }
    
    // HMMER approach: Log-odds scoring (HMM vs null model)
    // This is the statistically correct method for sequence alignment
    double null_score = 0.0;
    for (const auto& aa : aa_sequence) {
        int aa_idx = AminoAcidToIndex(aa[0]);
        if (aa_idx >= 0 && aa_idx < static_cast<int>(compo_match_.size())) {
            null_score += compo_match_[aa_idx];  // Background frequency (log-prob)
        } else {
            // Stop codons should never reach here (filtered by -INFINITY)
            // But handle gracefully if they do
            null_score += compo_match_[0];  // Use some default (e.g., 'A')
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
