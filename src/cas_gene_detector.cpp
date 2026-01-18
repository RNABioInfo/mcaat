#include "cas_gene_detector.h"
#include <algorithm>
#include <unordered_map>
#include <cmath>

CasGeneDetector::CasGeneDetector(SDBG& sdbg) : sdbg(sdbg), profile_(nullptr) {}
CasGeneDetector::CasGeneDetector(SDBG& sdbg, Profile* profile) : sdbg(sdbg), profile_(profile) {}

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
    int L = profile_->GetLength();
    col.Initialize(L);
    return col;
}

ViterbiColumn CasGeneDetector::ExtendViterbi(const ViterbiColumn& prev, char aa) {
    ViterbiColumn curr;
    if (!profile_) return curr;
    
    int L = profile_->GetLength();
    if (L <= 0) return curr;
    
    curr.M.assign(L + 1, -1e9);
    curr.I.assign(L + 1, -1e9);
    curr.D.assign(L + 1, -1e9);
    curr.seq_length = prev.seq_length + 1;
    curr.best_score = -1e9;  // Reset to find new best
    curr.best_hmm_pos = prev.best_hmm_pos;
    
    // Match states
    for (int j = 1; j <= L; ++j) {
        double emit_m = profile_->GetMatchEmission(j, aa);
        
        // For j=1, can enter from M[0] (begin state) with score 0
        // For j>1, normal transitions from previous states at j-1
        double from_m, from_i, from_d;
        
        if (j == 1) {
            // Entry from begin state - M[0] represents B->M1
            from_m = prev.M[0];  // prev.M[0] = 0 initially (begin)
            from_i = -1e9;       // Can't come from I at position 0
            from_d = -1e9;       // Can't come from D at position 0
        } else {
            from_m = prev.M[j-1] + profile_->GetTransition(j-1, j, 'M', 'M');
            from_i = prev.I[j-1] + profile_->GetTransition(j-1, j, 'I', 'M');
            from_d = prev.D[j-1] + profile_->GetTransition(j-1, j, 'D', 'M');
        }
        
        curr.M[j] = emit_m + std::max({from_m, from_i, from_d});
        
        // Track best score (for local alignment, any M state can be best)
        if (curr.M[j] > curr.best_score) {
            curr.best_score = curr.M[j];
            curr.best_hmm_pos = j;
        }
    }
    
    // Insert states (emit then stay or transition)
    for (int j = 1; j <= L; ++j) {
        double emit_i = profile_->GetInsertEmission(j, aa);
        double from_m = prev.M[j] + profile_->GetTransition(j, j, 'M', 'I');
        double from_i = prev.I[j] + profile_->GetTransition(j, j, 'I', 'I');
        curr.I[j] = emit_i + std::max(from_m, from_i);
    }
    
    // Delete states (no emission, computed from current row)
    for (int j = 2; j <= L; ++j) {
        double from_m = curr.M[j-1] + profile_->GetTransition(j-1, j, 'M', 'D');
        double from_d = curr.D[j-1] + profile_->GetTransition(j-1, j, 'D', 'D');
        curr.D[j] = std::max(from_m, from_d);
    }
    
    // Keep best score if current is worse (monotonic in local alignment)
    if (prev.best_score > curr.best_score) {
        curr.best_score = prev.best_score;
        curr.best_hmm_pos = prev.best_hmm_pos;
    }
    
    return curr;
}

std::vector<AminoAcidPathInfo> CasGeneDetector::BeamSearchAminoAcids(
    uint64_t start_node, int beam_width, int max_depth, int start_codon_offset) {
    
    std::vector<AminoAcidPathInfo> results;
    if (!sdbg.IsValidEdge(start_node)) return results;
    
    struct State {
        uint64_t node;
        std::string dna;
        std::string aa;
        ViterbiColumn vit;
        uint64_t path_hash;  // Rolling hash for cycle detection
        double score() const { return vit.best_score; }
    };
    
    // Pre-allocate
    std::vector<State> beam;
    beam.reserve(beam_width * 4);
    std::vector<State> next_beam;
    next_beam.reserve(beam_width * 4);
    
    // Initialize
    State init;
    init.node = start_node;
    init.dna = GetNodeSequence(start_node).substr(start_codon_offset);
    init.dna.reserve(max_depth + 50);
    init.aa.reserve(max_depth / 3 + 20);
    init.vit = InitializeViterbi();
    init.path_hash = start_node;
    
    // Score initial codons
    size_t init_codons = init.dna.length() / 3;
    for (size_t i = 0; i < init_codons; ++i) {
        char aa = CodonToAminoAcid(init.dna.substr(i * 3, 3));
        if (aa != '*') {
            init.aa += aa;
            if (profile_) init.vit = ExtendViterbi(init.vit, aa);
        }
    }
    
    beam.push_back(std::move(init));
    
    // Traverse - max depth is HMM length + 25% buffer (in bp = aa * 3)
    int effective_max_depth = profile_ ? static_cast<int>(profile_->GetLength() * 1.25 * 3) : max_depth;
    
    for (int depth = 0; depth < effective_max_depth && !beam.empty(); ++depth) {
        next_beam.clear();
        
        for (auto& s : beam) {
            uint64_t out[4];
            int outdeg = sdbg.OutgoingEdges(s.node, out);
            
            // Check termination: no outgoing edges OR HMM complete
            if (outdeg == 0 || (profile_ && s.vit.best_hmm_pos >= profile_->GetLength())) {
                // Terminal - save result
                AminoAcidPathInfo r;
                r.dna_sequence = std::move(s.dna);
                r.total_score = s.vit.best_score;
                r.hmm_position = s.vit.best_hmm_pos;
                r.is_complete = (profile_ && s.vit.best_hmm_pos >= profile_->GetLength());
                for (char c : s.aa) r.amino_acids.push_back(std::string(1, c));
                // No node path stored - DNA sequence is sufficient
                results.push_back(std::move(r));
                continue;
            }
            
            for (int i = 0; i < outdeg; ++i) {
                if (!sdbg.IsValidEdge(out[i])) continue;
                
                // Simple cycle check - just avoid immediate self-loops
                if (out[i] == s.node) continue;
                
                // Get the new nucleotide
                char new_nuc = GetNodeSequence(out[i]).back();
                
                // Check if this creates a new codon
                size_t new_dna_len = s.dna.length() + 1;
                size_t num_codons = new_dna_len / 3;
                size_t prev_codons = s.dna.length() / 3;
                
                State ns;
                ns.node = out[i];
                ns.path_hash = s.path_hash ^ (out[i] * 0x9e3779b97f4a7c15ULL);
                
                // Build DNA (append single char is cheap)
                ns.dna = s.dna;
                ns.dna += new_nuc;
                
                // Only process if new codon formed
                if (num_codons > prev_codons) {
                    std::string codon = ns.dna.substr(prev_codons * 3, 3);
                    char aa = CodonToAminoAcid(codon);
                    
                    ns.aa = s.aa;
                    ns.vit = s.vit;  // Copy Viterbi only when extending
                    
                    if (aa != '*') {
                        ns.aa += aa;
                        if (profile_) ns.vit = ExtendViterbi(ns.vit, aa);
                    }
                } else {
                    // No new codon - just copy references (shallow copy is fine)
                    ns.aa = s.aa;
                    ns.vit = s.vit;
                }
                
                next_beam.push_back(std::move(ns));
            }
        }
        
        // Prune beam
        if (next_beam.size() > static_cast<size_t>(beam_width)) {
            std::partial_sort(next_beam.begin(), next_beam.begin() + beam_width, next_beam.end(),
                [](const State& a, const State& b) { return a.score() > b.score(); });
            next_beam.resize(beam_width);
        }
        
        std::swap(beam, next_beam);
    }
    
    // Save remaining paths (reached max_depth)
    for (auto& s : beam) {
        AminoAcidPathInfo r;
        r.dna_sequence = std::move(s.dna);
        r.total_score = s.vit.best_score;
        r.hmm_position = s.vit.best_hmm_pos;
        r.is_complete = false;
        for (char c : s.aa) r.amino_acids.push_back(std::string(1, c));
        results.push_back(std::move(r));
    }
    
    // Sort results by score
    std::sort(results.begin(), results.end(), 
              [](const AminoAcidPathInfo& a, const AminoAcidPathInfo& b) {
                  return a.total_score > b.total_score;
              });
    
    return results;
}
