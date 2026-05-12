#include "cas/cas_gene_detector.h"
#include "cas/cas_workflow.h"  // For SearchDirection enum
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <cmath>

CasGeneDetector::CasGeneDetector(SDBG& sdbg) : sdbg(sdbg), profile_(nullptr) {}
CasGeneDetector::CasGeneDetector(SDBG& sdbg, const Profile* profile) : sdbg(sdbg), profile_(profile) {}

char CasGeneDetector::CodonToAminoAcid(const std::string& codon) {
    static const std::unordered_map<std::string, char> table = {
        {"TTT",'F'},{"TTC",'F'},{"TTA",'L'},{"TTG",'L'},
        {"TCT",'S'},{"TCC",'S'},{"TCA",'S'},{"TCG",'S'},
        {"TAT",'Y'},{"TAC",'Y'},{"TAA",'*'},{"TAG",'*'},
        {"TGT",'C'},{"TGC",'C'},{"TGA",'*'},{"TGG",'W'},
        {"CTT",'L'},{"CTC",'L'},{"CTA",'L'},{"CTG",'L'},
        {"CCT",'P'},{"CCC",'P'},{"CCA",'P'},{"CCG",'P'},
        {"CAT",'H'},{"CAC",'H'},{"CAA",'Q'},{"CAG",'Q'},
        {"CGT",'R'},{"CGC",'R'},{"CGA",'R'},{"CGG",'R'},
        {"ATT",'I'},{"ATC",'I'},{"ATA",'I'},{"ATG",'M'},
        {"ACT",'T'},{"ACC",'T'},{"ACA",'T'},{"ACG",'T'},
        {"AAT",'N'},{"AAC",'N'},{"AAA",'K'},{"AAG",'K'},
        {"AGT",'S'},{"AGC",'S'},{"AGA",'R'},{"AGG",'R'},
        {"GTT",'V'},{"GTC",'V'},{"GTA",'V'},{"GTG",'V'},
        {"GCT",'A'},{"GCC",'A'},{"GCA",'A'},{"GCG",'A'},
        {"GAT",'D'},{"GAC",'D'},{"GAA",'E'},{"GAG",'E'},
        {"GGT",'G'},{"GGC",'G'},{"GGA",'G'},{"GGG",'G'}
    };
    if (codon.length() != 3) return 'X';
    auto it = table.find(codon);
    return (it != table.end()) ? it->second : 'X';
}

std::string CasGeneDetector::GetNodeSequence(uint64_t node_id) {
    const uint32_t k = sdbg.k();
    std::vector<uint8_t> seq(k);
    sdbg.GetLabel(node_id, seq.data());
    std::string result;
    result.reserve(k);
    for (uint32_t i = 0; i < k; ++i) {
        switch (seq[i]) {
            case 1: result += 'A'; break;
            case 2: result += 'C'; break;
            case 3: result += 'G'; break;
            case 4: result += 'T'; break;
            default: result += 'N'; break;
        }
    }
    return result;
}

ViterbiColumn CasGeneDetector::InitializeViterbi() {
    ViterbiColumn col;
    if (!profile_) return col;
    col.Initialize(profile_->GetLength());
    return col;
}

ViterbiColumn CasGeneDetector::ExtendViterbiBanded(const ViterbiColumn& prev, char aa) {
    ViterbiColumn curr;
    if (!profile_) return curr;

    const int L = std::min(profile_->GetLength(), static_cast<int>(profile_->GetStates().size()));
    if (L <= 0) return curr;

    curr.M.assign(L + 1, -1e9);
    curr.I.assign(L + 1, -1e9);
    curr.D.assign(L + 1, -1e9);
    curr.E = -1e9;
    curr.seq_length = prev.seq_length + 1;
    curr.best_score = -1e9;
    curr.best_hmm_pos = prev.best_hmm_pos;
    
    // Band center follows best HMM position from previous column
    // For first few AAs, use wider band to allow entry at any position
    int center = (prev.seq_length < BAND_WIDTH) ? (L / 2) : (prev.best_hmm_pos + 1);
    curr.band_center = center;

    const double esc = profile_->IsLocal() ? 0.0 : -1e9;
    curr.M[0] = curr.I[0] = curr.D[0] = -1e9;

    // Banded computation: for early positions, use full range
    int k_min, k_max;
    if (prev.seq_length < BAND_WIDTH) {
        // Early: scan full profile to find best entry point
        k_min = 1;
        k_max = L;
    } else {
        // Later: use band around current best position
        k_min = std::max(1, center - BAND_WIDTH);
        k_max = std::min(L, center + BAND_WIDTH);
    }

    for (int k = k_min; k <= k_max; ++k) {
        double emission = profile_->GetMatchLogOdds(k, aa);

        // Match state
        if (prev.seq_length == 0) {
            curr.M[k] = prev.M[0] + emission;
        } else {
            double sc_m = prev.M[k-1] + profile_->GetTransition(k-1, k, 'M', 'M');
            sc_m = std::max(sc_m, prev.I[k-1] + profile_->GetTransition(k-1, k, 'I', 'M'));
            sc_m = std::max(sc_m, prev.D[k-1] + profile_->GetTransition(k-1, k, 'D', 'M'));
            curr.M[k] = sc_m + emission;
        }

        // Insert state
        double sc_i = prev.M[k] + profile_->GetTransition(k, k, 'M', 'I');
        sc_i = std::max(sc_i, prev.I[k] + profile_->GetTransition(k, k, 'I', 'I'));
        curr.I[k] = sc_i + profile_->GetInsertLogOdds(k, aa);

        // E state update
        curr.E = std::max(curr.E, curr.M[k] + esc);

        if (curr.M[k] > curr.best_score) {
            curr.best_score = curr.M[k];
            curr.best_hmm_pos = k;
        }
    }

    // Delete states: process separately (no emission, chain within current column)
    // Start from k_min-1 to allow D state to propagate into band
    int d_start = std::max(1, k_min - 1);
    for (int k = d_start; k <= k_max; ++k) {
        double sc_d = curr.M[k-1] + profile_->GetTransition(k-1, k, 'M', 'D');
        sc_d = std::max(sc_d, curr.D[k-1] + profile_->GetTransition(k-1, k, 'D', 'D'));
        curr.D[k] = sc_d;
    }

    curr.best_score = std::max(curr.best_score, curr.E);
    return curr;
}

std::vector<AminoAcidPathInfo> CasGeneDetector::BeamSearchAminoAcids(
    uint64_t start_node, int beam_width, int max_depth, int start_codon_offset,
    SearchDirection direction, double min_normalized_score) {
    
    std::vector<AminoAcidPathInfo> results;
    if (!sdbg.IsValidEdge(start_node)) return results;
    
    // NOTE: Gene traversal is ALWAYS forward (START→STOP) using OutgoingEdges!
    // The 'direction' parameter only affects how we FIND candidates, not how we traverse genes.
    // Genes transcribe 5'→3' regardless of their position relative to the repeat.
    (void)direction;  // Direction is irrelevant for gene body traversal
    
    // Pre-compute max emission score for early termination bound
    double max_emission_per_pos = profile_ ? profile_->GetMaxEmissionScore() : 3.0;
    int hmm_len = profile_ ? profile_->GetLength() : (max_depth / 3);
    
    struct State {
        uint64_t node;
        std::string dna;
        std::string aa;
        ViterbiColumn vit;
        std::vector<uint64_t> path;
        double score() const { return vit.best_score; }
    };
    
    std::vector<State> beam, next_beam;
    beam.reserve(beam_width * 4);
    next_beam.reserve(beam_width * 4);
    
    // Initialize
    State init;
    init.node = start_node;
    init.dna = GetNodeSequence(start_node).substr(start_codon_offset);
    init.vit = InitializeViterbi();
    init.path.push_back(start_node);

    // Score initial codons from START codon forward
    // Always: START codon → gene body → STOP codon (5'→3' direction)
    for (size_t i = 0; i + 2 < init.dna.length(); i += 3) {
        char aa = CodonToAminoAcid(init.dna.substr(i, 3));
        if (aa == '*') break;  // Stop at STOP codon
        init.aa += aa;
        if (profile_) init.vit = ExtendViterbiBanded(init.vit, aa);
    }
    beam.push_back(std::move(init));
    
    // Max traversal depth (hmm_len already computed at top for early termination)
    int max_aa = static_cast<int>(std::ceil(hmm_len * 1.25));
    int effective_max_depth = max_aa * 3;

    auto save_result = [&](State& s, bool complete) {
        if (s.aa.length() < 30) return;  // Too short
        int L = profile_ ? profile_->GetLength() : 1;
        AminoAcidPathInfo r;
        
        // No reversal needed - we always collect forward (START→STOP)
        r.dna_sequence = std::move(s.dna);
        r.total_score = s.vit.best_score / std::log(2.0);
        r.normalized_score = (L > 0) ? (r.total_score / L) : 0.0;
        r.hmm_position = s.vit.best_hmm_pos;
        r.is_complete = complete;
        r.node_path = std::move(s.path);
        for (char c : s.aa) r.amino_acids.push_back(std::string(1, c));
        results.push_back(std::move(r));
    };
    
    for (int depth = 0; depth < effective_max_depth && !beam.empty(); ++depth) {
        next_beam.clear();

        for (auto& s : beam) {
            uint64_t edges[4];
            int edge_count;
            
            // ALWAYS use OutgoingEdges - genes transcribe 5'→3' (START→STOP)
            edge_count = sdbg.OutgoingEdges(s.node, edges);

            // Dead end or HMM complete
            if (edge_count == 0 || (profile_ && s.vit.best_hmm_pos >= profile_->GetLength())) {
                save_result(s, true);
                continue;
            }
            
            for (int i = 0; i < edge_count; ++i) {
                if (!sdbg.IsValidEdge(edges[i]) || edges[i] == s.node) continue;
                
                // For outgoing edge, we want the LAST nucleotide of the next node
                char new_nuc = GetNodeSequence(edges[i]).back();
                
                size_t new_len = s.dna.length() + 1;
                size_t num_codons = new_len / 3;
                size_t prev_codons = s.dna.length() / 3;

                if (profile_ && num_codons > static_cast<size_t>(max_aa)) continue;
                
                State ns;
                ns.node = edges[i];
                ns.dna = s.dna + new_nuc;
                ns.path = s.path;
                ns.path.push_back(edges[i]);
                ns.aa = s.aa;
                ns.vit = s.vit;
                
                // New codon formed
                if (num_codons > prev_codons) {
                    std::string codon = ns.dna.substr(prev_codons * 3, 3);
                    char aa = CodonToAminoAcid(codon);

                    // STOP codon terminates the gene (always - both directions)
                    if (aa == '*') {
                        save_result(ns, true);
                        continue;
                    }
                    
                    ns.aa += aa;
                    if (profile_) ns.vit = ExtendViterbiBanded(ns.vit, aa);
                    
                    // EARLY TERMINATION: Check if this path can possibly reach threshold
                    // Only check every 10 amino acids to reduce overhead
                    // DISABLED FOR DEBUG - uncomment to enable
                    /*
                    if (profile_ && ns.aa.size() % 10 == 0 && ns.aa.size() >= 20) {
                        double current_score = ns.vit.best_score / std::log(2.0);  // Convert to bits
                        int remaining_pos = hmm_len - ns.vit.best_hmm_pos;
                        if (remaining_pos > 0) {
                            double best_possible = current_score + remaining_pos * max_emission_per_pos;
                            double best_normalized = best_possible / hmm_len;
                            if (best_normalized < min_normalized_score) {
                                // Mathematically impossible to reach threshold - skip this path
                                continue;
                            }
                        }
                    }
                    */
                }

                next_beam.push_back(std::move(ns));
            }
        }
        
        // Beam pruning
        if (next_beam.size() > static_cast<size_t>(beam_width)) {
            std::partial_sort(next_beam.begin(), next_beam.begin() + beam_width, next_beam.end(),
                [](const State& a, const State& b) { return a.score() > b.score(); });
            next_beam.resize(beam_width);
        }
        
        std::swap(beam, next_beam);
    }
    
    // Save remaining paths
    for (auto& s : beam) save_result(s, false);
    
    std::sort(results.begin(), results.end(), 
        [](const AminoAcidPathInfo& a, const AminoAcidPathInfo& b) {
            return a.total_score > b.total_score;
        });
    
    return results;
}

// ============================================================================
// Multi-Profile Viterbi Helpers
// ============================================================================

ViterbiColumn CasGeneDetector::InitializeViterbi(int hmm_length) {
    ViterbiColumn col;
    col.Initialize(hmm_length);
    return col;
}

ViterbiColumn CasGeneDetector::ExtendViterbiBanded(const ViterbiColumn& prev, char aa, const Profile* profile) {
    ViterbiColumn curr;
    if (!profile) return curr;

    const int L = std::min(profile->GetLength(), static_cast<int>(profile->GetStates().size()));
    if (L <= 0) return curr;

    curr.M.assign(L + 1, -1e9);
    curr.I.assign(L + 1, -1e9);
    curr.D.assign(L + 1, -1e9);
    curr.E = -1e9;
    curr.seq_length = prev.seq_length + 1;
    curr.best_score = -1e9;
    curr.best_hmm_pos = prev.best_hmm_pos;
    
    int center = (prev.seq_length < BAND_WIDTH) ? (L / 2) : (prev.best_hmm_pos + 1);
    curr.band_center = center;

    const double esc = profile->IsLocal() ? 0.0 : -1e9;
    curr.M[0] = curr.I[0] = curr.D[0] = -1e9;

    int k_min, k_max;
    if (prev.seq_length < BAND_WIDTH) {
        k_min = 1;
        k_max = L;
    } else {
        k_min = std::max(1, center - BAND_WIDTH);
        k_max = std::min(L, center + BAND_WIDTH);
    }

    for (int k = k_min; k <= k_max; ++k) {
        double emission = profile->GetMatchLogOdds(k, aa);

        if (prev.seq_length == 0) {
            curr.M[k] = prev.M[0] + emission;
        } else {
            double sc_m = prev.M[k-1] + profile->GetTransition(k-1, k, 'M', 'M');
            sc_m = std::max(sc_m, prev.I[k-1] + profile->GetTransition(k-1, k, 'I', 'M'));
            sc_m = std::max(sc_m, prev.D[k-1] + profile->GetTransition(k-1, k, 'D', 'M'));
            curr.M[k] = sc_m + emission;
        }

        double sc_i = prev.M[k] + profile->GetTransition(k, k, 'M', 'I');
        sc_i = std::max(sc_i, prev.I[k] + profile->GetTransition(k, k, 'I', 'I'));
        curr.I[k] = sc_i + profile->GetInsertLogOdds(k, aa);

        curr.E = std::max(curr.E, curr.M[k] + esc);

        if (curr.M[k] > curr.best_score) {
            curr.best_score = curr.M[k];
            curr.best_hmm_pos = k;
        }
    }

    int d_start = std::max(1, k_min - 1);
    for (int k = d_start; k <= k_max; ++k) {
        double sc_d = curr.M[k-1] + profile->GetTransition(k-1, k, 'M', 'D');
        sc_d = std::max(sc_d, curr.D[k-1] + profile->GetTransition(k-1, k, 'D', 'D'));
        curr.D[k] = sc_d;
    }

    curr.best_score = std::max(curr.best_score, curr.E);
    return curr;
}

// ============================================================================
// Multi-Profile Beam Search: ONE traversal, ALL profiles scored simultaneously
// ============================================================================

std::vector<MultiProfileResult> CasGeneDetector::BeamSearchMultiProfile(
    uint64_t start_node, int beam_width, int max_depth, int start_codon_offset,
    const std::vector<std::pair<const Profile*, size_t>>& profiles,
    double min_normalized_score) {
    
    std::vector<MultiProfileResult> results;
    if (!sdbg.IsValidEdge(start_node) || profiles.empty()) return results;
    
    const size_t num_profiles = profiles.size();
    
    // Find max HMM length for max_aa calculation
    int max_hmm_len = 0;
    for (const auto& [profile, idx] : profiles) {
        if (profile) max_hmm_len = std::max(max_hmm_len, profile->GetLength());
    }
    int max_aa = static_cast<int>(std::ceil(max_hmm_len * 1.25));
    int effective_max_depth = max_aa * 3;
    
    // State now carries Viterbi columns for ALL profiles
    struct MultiState {
        uint64_t node;
        std::string dna;
        std::string aa;
        std::vector<ViterbiColumn> viterbi;  // One per profile
        std::vector<uint64_t> path;
        
        // Max score across all profiles (for beam pruning)
        double max_score() const {
            double m = -1e9;
            for (const auto& v : viterbi) {
                m = std::max(m, v.best_score);
            }
            return m;
        }
        
        // Score for a specific profile
        double score(size_t p) const {
            return (p < viterbi.size()) ? viterbi[p].best_score : -1e9;
        }
    };
    
    std::vector<MultiState> beam, next_beam;
    beam.reserve(beam_width * 4);
    next_beam.reserve(beam_width * 4);
    
    // Initialize
    MultiState init;
    init.node = start_node;
    init.dna = GetNodeSequence(start_node).substr(start_codon_offset);
    init.path.push_back(start_node);
    
    // Initialize Viterbi for each profile
    init.viterbi.resize(num_profiles);
    for (size_t p = 0; p < num_profiles; ++p) {
        const Profile* profile = profiles[p].first;
        if (profile) {
            init.viterbi[p] = InitializeViterbi(profile->GetLength());
        }
    }
    
    // Score initial codons from START codon forward
    for (size_t i = 0; i + 2 < init.dna.length(); i += 3) {
        char aa = CodonToAminoAcid(init.dna.substr(i, 3));
        if (aa == '*') break;
        init.aa += aa;
        for (size_t p = 0; p < num_profiles; ++p) {
            const Profile* profile = profiles[p].first;
            if (profile) {
                init.viterbi[p] = ExtendViterbiBanded(init.viterbi[p], aa, profile);
            }
        }
    }
    beam.push_back(std::move(init));
    
    auto save_result = [&](MultiState& s, bool complete) {
        if (s.aa.length() < 30) return;
        
        MultiProfileResult r;
        r.dna_sequence = std::move(s.dna);
        r.amino_acids = std::move(s.aa);
        r.node_path = std::move(s.path);
        r.is_complete = complete;
        
        // Compute scores for each profile
        r.profile_scores.resize(num_profiles);
        for (size_t p = 0; p < num_profiles; ++p) {
            const Profile* profile = profiles[p].first;
            int L = profile ? profile->GetLength() : 1;
            double total = s.viterbi[p].best_score / std::log(2.0);
            double normalized = (L > 0) ? (total / L) : 0.0;
            int hmm_pos = s.viterbi[p].best_hmm_pos;
            r.profile_scores[p] = std::make_tuple(total, normalized, hmm_pos);
        }
        
        results.push_back(std::move(r));
    };
    
    for (int depth = 0; depth < effective_max_depth && !beam.empty(); ++depth) {
        next_beam.clear();
        
        for (auto& s : beam) {
            uint64_t edges[4];
            int edge_count = sdbg.OutgoingEdges(s.node, edges);
            
            // Dead end or all profiles complete
            if (edge_count == 0) {
                save_result(s, true);
                continue;
            }
            
            for (int i = 0; i < edge_count; ++i) {
                if (!sdbg.IsValidEdge(edges[i]) || edges[i] == s.node) continue;
                
                char new_nuc = GetNodeSequence(edges[i]).back();
                
                size_t new_len = s.dna.length() + 1;
                size_t num_codons = new_len / 3;
                size_t prev_codons = s.dna.length() / 3;
                
                if (num_codons > static_cast<size_t>(max_aa)) continue;
                
                MultiState ns;
                ns.node = edges[i];
                ns.dna = s.dna + new_nuc;
                ns.path = s.path;
                ns.path.push_back(edges[i]);
                ns.aa = s.aa;
                ns.viterbi = s.viterbi;  // Copy all Viterbi states
                
                // New codon formed
                if (num_codons > prev_codons) {
                    std::string codon = ns.dna.substr(prev_codons * 3, 3);
                    char aa = CodonToAminoAcid(codon);
                    
                    if (aa == '*') {
                        save_result(ns, true);
                        continue;
                    }
                    
                    ns.aa += aa;
                    
                    // Score this AA against ALL profiles
                    for (size_t p = 0; p < num_profiles; ++p) {
                        const Profile* profile = profiles[p].first;
                        if (profile) {
                            ns.viterbi[p] = ExtendViterbiBanded(ns.viterbi[p], aa, profile);
                        }
                    }
                }
                
                next_beam.push_back(std::move(ns));
            }
        }
        
        // INFORMATION-THEORETIC BEAM PRUNING
        // Keep state if it's in top-K for ANY profile
        if (next_beam.size() > static_cast<size_t>(beam_width)) {
            // For each profile, mark top-K states
            std::vector<bool> keep(next_beam.size(), false);
            
            for (size_t p = 0; p < num_profiles; ++p) {
                // Get indices sorted by this profile's score
                std::vector<size_t> indices(next_beam.size());
                std::iota(indices.begin(), indices.end(), 0);
                std::partial_sort(indices.begin(), 
                    indices.begin() + std::min(static_cast<size_t>(beam_width), indices.size()),
                    indices.end(),
                    [&](size_t a, size_t b) {
                        return next_beam[a].score(p) > next_beam[b].score(p);
                    });
                
                // Mark top-K for this profile
                for (size_t k = 0; k < std::min(static_cast<size_t>(beam_width), indices.size()); ++k) {
                    keep[indices[k]] = true;
                }
            }
            
            // Filter to kept states
            std::vector<MultiState> filtered;
            filtered.reserve(beam_width * num_profiles);
            for (size_t i = 0; i < next_beam.size(); ++i) {
                if (keep[i]) {
                    filtered.push_back(std::move(next_beam[i]));
                }
            }
            next_beam = std::move(filtered);
            
            // If still too large, fall back to max-score pruning
            if (next_beam.size() > static_cast<size_t>(beam_width * 3)) {
                std::partial_sort(next_beam.begin(),
                    next_beam.begin() + beam_width * 2,
                    next_beam.end(),
                    [](const MultiState& a, const MultiState& b) {
                        return a.max_score() > b.max_score();
                    });
                next_beam.resize(beam_width * 2);
            }
        }
        
        std::swap(beam, next_beam);
    }
    
    // Save remaining
    for (auto& s : beam) save_result(s, false);
    
    return results;
}