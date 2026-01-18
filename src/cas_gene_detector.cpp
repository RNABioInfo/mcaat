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

    const int L = std::min(profile_->GetLength(), static_cast<int>(profile_->GetStates().size()));
    if (L <= 0) return curr;

    curr.M.assign(L + 1, -1e9);
    curr.I.assign(L + 1, -1e9);
    curr.D.assign(L + 1, -1e9);
    curr.E = -1e9;
    curr.seq_length = prev.seq_length + 1;
    curr.best_score = -1e9;
    curr.best_hmm_pos = prev.best_hmm_pos;

    // HMMER-style local esc: 0 for local, -inf for glocal
    const double esc = profile_->IsLocal() ? 0.0 : -1e9;

    // Set row 0
    curr.M[0] = curr.I[0] = curr.D[0] = -1e9;

    // Core recurrences (match/insert/delete), no B/N/J/C extras to keep it tight
    for (int k = 1; k <= L; ++k) {
        // Match
        double sc_m = prev.M[k-1] + profile_->GetTransition(k-1, k, 'M', 'M');
        sc_m = std::max(sc_m, prev.I[k-1] + profile_->GetTransition(k-1, k, 'I', 'M'));
        sc_m = std::max(sc_m, prev.D[k-1] + profile_->GetTransition(k-1, k, 'D', 'M'));
        curr.M[k] = sc_m + profile_->GetMatchLogOdds(k, aa);

        // Insert
        double sc_i = prev.M[k] + profile_->GetTransition(k, k, 'M', 'I');
        sc_i = std::max(sc_i, prev.I[k] + profile_->GetTransition(k, k, 'I', 'I'));
        curr.I[k] = sc_i + profile_->GetInsertLogOdds(k, aa);

        // Delete (depends on current row M/D at k-1)
        double sc_d = curr.M[k-1] + profile_->GetTransition(k-1, k, 'M', 'D');
        sc_d = std::max(sc_d, curr.D[k-1] + profile_->GetTransition(k-1, k, 'D', 'D'));
        curr.D[k] = sc_d;

        // E update (local): allow only from M as in HMMER local path
        curr.E = std::max(curr.E, curr.M[k] + esc);

        if (curr.M[k] > curr.best_score) {
            curr.best_score = curr.M[k];
            curr.best_hmm_pos = k;
        }
    }

    // Best score for this row: E (local) or best M if glocal
    curr.best_score = std::max(curr.best_score, curr.E);
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
        std::vector<uint64_t> path;  // node path
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
    init.path.push_back(start_node);
    
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
    
    // Traverse - max depth is HMM length + 25% buffer (aa -> bp)
    int max_aa = profile_ ? static_cast<int>(std::ceil(std::min(profile_->GetLength(), static_cast<int>(profile_->GetStates().size())) * 1.25)) : max_depth / 3;
    int effective_max_depth = max_aa * 3;
    
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
                // Convert log-odds (nats) to bits: bits = nats / ln(2)
                r.total_score = s.vit.best_score / std::log(2.0);
                r.hmm_position = s.vit.best_hmm_pos;
                r.is_complete = (profile_ && s.vit.E >= s.vit.best_score);
                r.node_path = std::move(s.path);
                for (char c : s.aa) r.amino_acids.push_back(std::string(1, c));
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

                // Stop if AA length exceeds profile length +25%
                if (profile_ && num_codons > static_cast<size_t>(max_aa)) continue;
                
                State ns;
                ns.node = out[i];
                ns.path_hash = s.path_hash ^ (out[i] * 0x9e3779b97f4a7c15ULL);
                
                // Build DNA (append single char is cheap)
                ns.dna = s.dna;
                ns.dna += new_nuc;

                // Propagate path and append next node
                ns.path = s.path;
                ns.path.push_back(out[i]);
                
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
        // Convert log-odds (nats) to bits
        r.total_score = s.vit.best_score / std::log(2.0);
        r.hmm_position = s.vit.best_hmm_pos;
        r.is_complete = false;
        r.node_path = std::move(s.path);
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
