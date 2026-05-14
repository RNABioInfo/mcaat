#include "cas/cas_workflow.h"
#include "cas/cas_gene_detector.h"
#include <algorithm>
#include <atomic>
#include <mutex>
#include <set>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <omp.h>

// Genes that appear in many types — skip candidate-type narrowing when detected.
// These are adaptation genes (cas1/2/3/4/6) and regulatory accessory genes.
static const std::set<std::string> MULTI_MODEL_FAMILIES = {
    "cas1", "cas2", "cas3", "cas4", "cas6", "cas6e", "cas6f",
    "rt", "wyl", "deddh", "casr", "primpol", "pd"
};

CasWorkflow::CasWorkflow(SDBG& sdbg, const std::string& profiles_dir, const std::string& rules_csv_path)
    : sdbg_(sdbg), profiles_dir_(profiles_dir) {
    std::cout << __func__ << "\n";

    // Ensure OpenMP uses full CPU by default unless user explicitly configured thread count.
    if (std::getenv("OMP_NUM_THREADS") == nullptr) {
        omp_set_dynamic(0);
        omp_set_num_threads(omp_get_num_procs());
    }

    profiles_ = HMMProfiles::LoadProfilesFromDirectory(profiles_dir);
    
    // Initialize lazy-loading storage
    loaded_profiles_.resize(profiles_.size());
    profile_loaded_.resize(profiles_.size(), false);
    
    // Load rules from CSV and build indices
    // Build family index FIRST so it's available for first-gene lookup
    BuildFamilyIndex();
    LoadRulesFromCSV(rules_csv_path);
    BuildFirstGeneIndex();
}

void CasWorkflow::ApplySensitivityMode(bool enabled) {
    sensitivity_mode_ = enabled;

    if (!enabled) {
        tuning_ = WorkflowTuning{};
        return;
    }

    const bool preserve_exploratory = params_.allow_exploratory;

    params_.FIRST_GENE_MIN_DIST = 0;
    params_.FIRST_GENE_MAX_DIST = 2500;
    params_.MAX_START_CANDIDATES = 60000;
    params_.BEAM_WIDTH = 100;
    params_.MIN_NORMALIZED_SCORE = 0.2;
    params_.allow_exploratory = false;

    tuning_.first_gene_search_max_bp = 10000;
    tuning_.bfs_max_candidates = 40000;
    tuning_.max_cassette_bp = 60000;
    tuning_.intergenic_min_bp = -99;
    tuning_.intergenic_max_bp = 500;
    tuning_.max_genes = 20;
    tuning_.first_gene_bins = 100;
    tuning_.profile_window_count = 100;
    tuning_.shallow_threshold = 0.05;
    tuning_.phase1_high_confidence = 0.95;
    tuning_.phase2_high_confidence = 0.99;
}
// Build a profile-specific candidate window by distance binning.
// Picks at most one START candidate per bin across the allowed range,
// up to max_per_profile total candidates.
static std::vector<const StartCodonCandidate*> BuildBinnedCandidateWindow(
    const std::vector<StartCodonCandidate>& sorted_candidates,
    int allowed_min_dist, int allowed_max_dist,
    int max_per_profile) {

    std::vector<const StartCodonCandidate*> selected;
    if (sorted_candidates.empty() || max_per_profile <= 0 || allowed_min_dist > allowed_max_dist) {
        return selected;
    }

    selected.reserve(max_per_profile);

    auto lower = std::lower_bound(sorted_candidates.begin(), sorted_candidates.end(), allowed_min_dist,
        [](const StartCodonCandidate& c, int d) { return c.distance < d; });
    auto upper = std::upper_bound(sorted_candidates.begin(), sorted_candidates.end(), allowed_max_dist,
        [](int d, const StartCodonCandidate& c) { return d < c.distance; });

    if (lower == upper) {
        return selected;
    }

    const int span = std::max(1, allowed_max_dist - allowed_min_dist + 1);
    const int bin_width = std::max(1, (span + max_per_profile - 1) / max_per_profile);
    std::vector<bool> bin_taken(max_per_profile, false);

    for (auto it = lower; it != upper; ++it) {
        int bin = (it->distance - allowed_min_dist) / bin_width;
        if (bin < 0) continue;
        if (bin >= max_per_profile) bin = max_per_profile - 1;

        if (!bin_taken[bin]) {
            selected.push_back(&(*it));
            bin_taken[bin] = true;
            if (static_cast<int>(selected.size()) >= max_per_profile) {
                break;
            }
        }
    }

    return selected;
}

void CasWorkflow::LoadRulesFromCSV(const std::string& rules_csv_path) {
    std::cout << __func__ << "\n";
    std::ifstream file(rules_csv_path);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open rules file: " << rules_csv_path << "\n";
        return;
    }
    
    std::string line;
    int line_num = 0;
    int skipped = 0;
    int parsed = 0;
    
    while (std::getline(file, line)) {
        ++line_num;
        
        // Skip comments and header
        if (line.empty() || line[0] == '#') {
            skipped++;
            continue;
        }
        // Skip header line
        if (line.find("type_class") != std::string::npos && line.find("mandatory") != std::string::npos) {
            skipped++;
            continue;
        }
        
        CasTypeRule rule;
        
        // CSV format: type_class,"mandatory_profiles","accessory_profiles",min_mandatory,min_genes,max_genes
        // Find first comma (after type_class)
        size_t pos = 0;
        size_t first_comma = line.find(',');
        if (first_comma == std::string::npos) continue;
        
        rule.type_class = line.substr(0, first_comma);
        pos = first_comma + 1;
        
        // Helper lambda to parse a field with semicolon-separated gene names (NO QUOTES)
        // Gene names with | represent alternatives (any one satisfies)
        auto parseField = [&](std::vector<std::string>& genes) {
            size_t next_comma = line.find(',', pos);
            size_t end_pos = (next_comma != std::string::npos) ? next_comma : line.size();
            
            std::string field = line.substr(pos, end_pos - pos);
            if (!field.empty()) {
                std::stringstream ss(field);
                std::string gene;
                while (std::getline(ss, gene, ';')) {
                    if (!gene.empty()) {
                        genes.push_back(gene);
                    }
                }
            }
            pos = (next_comma != std::string::npos) ? next_comma + 1 : line.size();
        };
        
        // Parse mandatory_genes (type-defining mandatory genes)
        parseField(rule.mandatory_genes);
        
        // Parse accessory_genes (optional genes including adaptation)
        parseField(rule.accessory_genes);
        
        // Parse min_mandatory, min_genes, and max_genes
        std::string remaining = line.substr(pos);
        std::stringstream rs(remaining);
        std::string val;
        int idx = 0;
        while (std::getline(rs, val, ',')) {
            try {
                int v = std::stoi(val);
                if (idx == 0) rule.min_mandatory = v;
                else if (idx == 1) rule.min_genes = v;
                else if (idx == 2) rule.max_genes = v;
                else if (idx == 3) rule.db_count  = std::max(1, v);
                ++idx;
            } catch (...) {}
        }
        
        if (!rule.mandatory_genes.empty() || !rule.accessory_genes.empty()) {
            parsed++;
            
            // Extract gene families from profile names for family-based matching
            // Each entry may contain alternatives separated by |
            auto extractFamilies = [](const std::string& gene_entry) -> std::set<std::string> {
                std::set<std::string> families;
                if (gene_entry.find('|') != std::string::npos) {
                    std::stringstream ss(gene_entry);
                    std::string alt;
                    while (std::getline(ss, alt, '|')) {
                        if (!alt.empty()) {
                            families.insert(ExtractGeneFamily(alt));
                        }
                    }
                } else {
                    families.insert(ExtractGeneFamily(gene_entry));
                }
                return families;
            };
            
            for (const auto& gene_entry : rule.mandatory_genes) {
                rule.mandatory_families.push_back(extractFamilies(gene_entry));
            }
            for (const auto& gene_entry : rule.accessory_genes) {
                rule.accessory_families.push_back(extractFamilies(gene_entry));
            }
            type_rules_.push_back(rule);
        }
    }
    
}

void CasWorkflow::BuildFirstGeneIndex() {
    std::cout << __func__ << "\n";
    std::set<size_t> first_gene_set;
    
    // Helper lambda to process gene entries (handles | alternatives)
    auto processGeneEntry = [&](const std::string& gene_entry, const std::string& type_class) {
        std::vector<std::string> alternatives;
        if (gene_entry.find('|') != std::string::npos) {
            std::stringstream ss(gene_entry);
            std::string alt;
            while (std::getline(ss, alt, '|')) {
                if (!alt.empty()) {
                    alternatives.push_back(alt);
                }
            }
        } else {
            alternatives.push_back(gene_entry);
        }
        
        for (const auto& gene_name : alternatives) {
            // Extract FAMILY from the gene name (e.g., "cas3_I" -> "cas3")
            std::string family = ExtractGeneFamily(gene_name);
            auto it = family_to_profiles_.find(family);
            
            if (it != family_to_profiles_.end()) {
                for (size_t idx : it->second) {
                    first_gene_set.insert(idx);
                    // Track which types this gene can indicate (avoid duplicates)
                    auto& types = first_gene_to_types_[idx];
                    if (std::find(types.begin(), types.end(), type_class) == types.end()) {
                        types.push_back(type_class);
                    }
                }
            }
        }
    };
    
    for (const auto& rule : type_rules_) {
        // ALL mandatory genes can start a cassette (type-defining from XML)
        for (const auto& gene_entry : rule.mandatory_genes) {
            processGeneEntry(gene_entry, rule.type_class);
        }
        
        // ALL accessory genes can start a cassette (accessory from XML)
        // Includes adaptation genes (cas1/2/4) and other optional genes (cas6, cas3, etc.)
        for (const auto& gene_entry : rule.accessory_genes) {
            processGeneEntry(gene_entry, rule.type_class);
        }
    }
    
    first_gene_profile_indices_.assign(first_gene_set.begin(), first_gene_set.end());
    
    if (first_gene_profile_indices_.empty()) {
        std::cerr << "ERROR: No first-gene profiles loaded! Detection will fail!\n";
    }
    
    // Pre-load first-gene profiles for fast Phase 1 detection
    for (size_t idx : first_gene_profile_indices_) {
        if (loaded_profiles_[idx].LoadFromFile(profiles_dir_ + "/" + profiles_[idx].filename)) {
            profile_loaded_[idx] = true;
        }
    }
    
    // PRE-LOAD ALL PROFILES to avoid lock contention in parallel loops
    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < profiles_.size(); ++idx) {
        if (!profile_loaded_[idx]) {
            #pragma omp critical(profile_load)
            {
                if (!profile_loaded_[idx]) {
                    loaded_profiles_[idx].LoadFromFile(profiles_dir_ + "/" + profiles_[idx].filename);
                    profile_loaded_[idx] = true;
                }
            }
        }
    }
}

void CasWorkflow::BuildFamilyIndex() {
    family_to_profiles_.clear();
    for (size_t i = 0; i < profiles_.size(); ++i) {
        std::string family = ExtractGeneFamily(profiles_[i].filename);
        family_to_profiles_[family].push_back(i);
    }
}

const Profile* CasWorkflow::GetProfile(size_t profile_index) const {
    if (profile_index >= loaded_profiles_.size()) {
        return nullptr;
    }
    
    // Check if already loaded (fast path, no lock needed for read)
    if (profile_loaded_[profile_index]) {
        const Profile* p = &loaded_profiles_[profile_index];
        return (p->GetLength() > 0) ? p : nullptr;
    }
    
    // Lazy load with mutex protection
    #pragma omp critical(profile_load)
    {
        // Double-check inside critical section
        if (!profile_loaded_[profile_index]) {
            loaded_profiles_[profile_index].LoadFromFile(
                profiles_dir_ + "/" + profiles_[profile_index].filename);
            profile_loaded_[profile_index] = true;
        }
    }
    
    const Profile* p = &loaded_profiles_[profile_index];
    return (p->GetLength() > 0) ? p : nullptr;
}

std::string CasWorkflow::GetNodeSequence(uint64_t node_id) {
    const uint32_t k = sdbg_.k();
    std::vector<uint8_t> seq(k);
    sdbg_.GetLabel(node_id, seq.data());
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

int CasWorkflow::GetStartCodonOffset(uint64_t node_id) {
    std::string seq = GetNodeSequence(node_id);
    if (seq.length() < 3) return -1;
    for (size_t i = 0; i + 2 < seq.length(); ++i) {
        std::string codon = seq.substr(i, 3);
        if (codon == "ATG" || codon == "GTG" || codon == "TTG") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool CasWorkflow::HasStartCodon(uint64_t node_id) {
    return GetStartCodonOffset(node_id) >= 0;
}

int CasWorkflow::GetStopCodonOffset(uint64_t node_id) {
    std::string seq = GetNodeSequence(node_id);
    if (seq.length() < 3) return -1;
    for (size_t i = 0; i + 2 < seq.length(); ++i) {
        std::string codon = seq.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool CasWorkflow::HasStopCodon(uint64_t node_id) {
    return GetStopCodonOffset(node_id) >= 0;
}

// Estimate minimum ORF length from a start codon position (DOWNSTREAM direction)
// Returns length in bp until first stop codon (or max_search if no stop found)
// This is a quick filter - doesn't need to be perfect, just catch obviously short ORFs
int CasWorkflow::EstimateORFLength(uint64_t start_node, int start_offset, int max_search) {
    const uint32_t k = sdbg_.k();
    std::string current_seq = GetNodeSequence(start_node).substr(start_offset);
    int bp_count = current_seq.length();
    
    // Check initial sequence for stop codons
    for (size_t i = 0; i + 2 < current_seq.length(); i += 3) {
        std::string codon = current_seq.substr(i, 3);
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return i;  // Found stop codon
        }
    }
    
    // If initial sequence is long enough, return early
    if (bp_count >= max_search) return bp_count;
    
    // Simple BFS to extend and find stop codon
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, std::string>> q;
    q.push({start_node, current_seq});
    visited.insert(start_node);
    
    while (!q.empty() && bp_count < max_search) {
        auto [node, seq] = q.front();
        q.pop();
        
        uint64_t out[4];
        int outdeg = sdbg_.OutgoingEdges(node, out);
        
        // Follow first valid edge (greedy approximation)
        for (int i = 0; i < outdeg; ++i) {
            if (sdbg_.IsValidEdge(out[i]) && !visited.count(out[i])
                && (cycle_nodes_.empty() || !cycle_nodes_[out[i]])) {
                visited.insert(out[i]);
                char next_base = GetNodeSequence(out[i]).back();
                std::string extended = seq + next_base;
                bp_count++;
                
                // Check if we completed a codon
                if (extended.length() >= 3 && extended.length() % 3 == 0) {
                    std::string last_codon = extended.substr(extended.length() - 3, 3);
                    if (last_codon == "TAA" || last_codon == "TAG" || last_codon == "TGA") {
                        return bp_count;  // Found stop codon
                    }
                }
                
                if (bp_count >= max_search) return bp_count;
                
                q.push({out[i], extended.length() > 100 ? extended.substr(extended.length() - 100) : extended});
                break;  // Only follow one path (greedy)
            }
        }
    }
    
    return bp_count;  // No stop found within search range
}

// Estimate minimum ORF length from a START codon position (UPSTREAM direction)
// Searches BACKWARDS using IncomingEdges to find STOP codon
// Returns length in bp until first stop codon (or max_search if no stop found)
int CasWorkflow::EstimateORFLengthReverse(uint64_t start_node, int start_offset, int max_search) {
    const uint32_t k = sdbg_.k();
    std::string node_seq = GetNodeSequence(start_node);
    
    // Start with sequence before start codon (reverse direction)
    std::string current_seq;
    if (start_offset > 0) {
        current_seq = node_seq.substr(0, start_offset);
        // Reverse it since we're going backwards
        std::reverse(current_seq.begin(), current_seq.end());
    }
    
    int bp_count = current_seq.length();
    
    // Check initial sequence for stop codons (reading backwards)
    for (size_t i = 0; i + 2 < current_seq.length(); i += 3) {
        std::string codon = current_seq.substr(i, 3);
        // Reverse the codon to check (we're reading backwards)
        std::reverse(codon.begin(), codon.end());
        if (codon == "TAA" || codon == "TAG" || codon == "TGA") {
            return i;  // Found stop codon
        }
    }
    
    if (bp_count >= max_search) return bp_count;
    
    // BFS backwards using IncomingEdges
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, std::string>> q;
    q.push({start_node, current_seq});
    visited.insert(start_node);
    
    while (!q.empty() && bp_count < max_search) {
        auto [node, seq] = q.front();
        q.pop();
        
        uint64_t in[4];
        int indeg = sdbg_.IncomingEdges(node, in);
        
        // Follow first valid edge (greedy approximation)
        for (int i = 0; i < indeg; ++i) {
            if (sdbg_.IsValidEdge(in[i]) && !visited.count(in[i])
                && (cycle_nodes_.empty() || !cycle_nodes_[in[i]])) {
                visited.insert(in[i]);
                char prev_base = GetNodeSequence(in[i])[0];  // First nucleotide of incoming node
                std::string extended = seq + prev_base;
                bp_count++;
                
                // Check if we completed a codon
                if (extended.length() >= 3 && extended.length() % 3 == 0) {
                    std::string last_codon = extended.substr(extended.length() - 3, 3);
                    // Reverse to check (we collected backwards)
                    std::reverse(last_codon.begin(), last_codon.end());
                    if (last_codon == "TAA" || last_codon == "TAG" || last_codon == "TGA") {
                        return bp_count;  // Found stop codon
                    }
                }
                
                if (bp_count >= max_search) return bp_count;
                
                q.push({in[i], extended.length() > 100 ? extended.substr(extended.length() - 100) : extended});
                break;  // Only follow one path (greedy)
            }
        }
    }
    
    return bp_count;  // No stop found within search range
}

// Check if a profile's length range is compatible with available space
// CRITICAL GEOMETRY:
//   UPSTREAM:   [START]--gene-->[STOP] <--D-- [REPEAT]
//   DOWNSTREAM: [REPEAT] --D--> [START]--gene-->[STOP]
// 
// min_bp/max_bp already encode the -15%/+25% tolerance from LENG (see buckets.h).
// UPSTREAM: Gene extends forward, available space = D
//   Constraint: D ∈ [min_bp + 100, max_bp + 100]
// DOWNSTREAM: Gene extends away, no constraint
bool CasWorkflow::IsProfileCompatibleWithDistance(
    const HMMProfiles::ProfileSize& profile, int distance_from_repeat,
    SearchDirection direction) {
    
    if (direction == SearchDirection::DOWNSTREAM) {
        return true;  // Gene extends away from repeat
    }
    
    // UPSTREAM: use min_bp/max_bp directly (already ±15%/+25% of LENG)
    // +100 accounts for intergenic gap between gene end and repeat
    double min_valid = profile.min_bp + 100.0;
    double max_valid = profile.max_bp + 100.0;
    
    if (distance_from_repeat < min_valid || distance_from_repeat > max_valid) {
        return false;
    }
    
    return true;
}

std::vector<StartCodonCandidate> CasWorkflow::FindStartCodonCandidates(
    uint64_t repeat_node, int min_dist, int max_dist, int max_candidates,
    SearchDirection direction) {
    std::vector<StartCodonCandidate> candidates;
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    
    q.push({repeat_node, 0});
    visited.insert(repeat_node);
    
    while (!q.empty() && static_cast<int>(candidates.size()) < max_candidates) {
        auto [node, dist] = q.front();
        q.pop();
        
        if (dist <= max_dist) {
            // Both directions now search for START codons
            int offset = GetStartCodonOffset(node);
            
            if (offset >= 0) {
                int bp_dist = dist + offset;
                if (bp_dist >= min_dist && bp_dist <= max_dist) {
                    // Lightweight ORF check: reject candidates that hit a stop codon
                    // within the first 150 bp (~50 aa) — minimum sensible Cas gene size.
                    // EstimateORFLength max_search=300 keeps this cheap (greedy, no full BFS).
                    constexpr int MIN_ORF_BP = 150;
                    int orf_len = (direction == SearchDirection::DOWNSTREAM)
                        ? EstimateORFLength(node, offset, 300)
                        : EstimateORFLengthReverse(node, offset, 300);
                    if (orf_len >= MIN_ORF_BP) {
                        candidates.push_back({node, bp_dist, offset, orf_len});
                    }
                }
            }
        }
        
        if (dist < max_dist) {
            uint64_t edges[4];
            int edge_count;
            
            // Direction-aware edge selection
            if (direction == SearchDirection::DOWNSTREAM) {
                edge_count = sdbg_.OutgoingEdges(node, edges);
            } else {
                edge_count = sdbg_.IncomingEdges(node, edges);
            }
            
            for (int i = 0; i < edge_count; ++i) {
                if (sdbg_.IsValidEdge(edges[i]) && !visited.count(edges[i])
                    && (cycle_nodes_.empty() || !cycle_nodes_[edges[i]])) {
                    q.push({edges[i], dist + 1});
                    visited.insert(edges[i]);
                }
            }
        }
    }
    return candidates;
}

std::vector<StartCodonCandidate> CasWorkflow::FindStartCodonCandidatesFromNode(
    uint64_t origin_node, int min_dist, int max_dist, int max_candidates,
    SearchDirection direction) {
    std::vector<StartCodonCandidate> candidates;
    std::set<uint64_t> visited;
    std::queue<std::pair<uint64_t, int>> q;
    
    q.push({origin_node, 0});
    visited.insert(origin_node);
    
    while (!q.empty() && static_cast<int>(candidates.size()) < max_candidates) {
        auto [node, dist] = q.front();
        q.pop();
        
        if (dist >= min_dist && dist <= max_dist) {
            // Both directions search for START codons
            int offset = GetStartCodonOffset(node);
            
            if (offset >= 0) {
                constexpr int MIN_ORF_BP = 120;
                int orf_len = (direction == SearchDirection::DOWNSTREAM)
                    ? EstimateORFLength(node, offset, 300)
                    : EstimateORFLengthReverse(node, offset, 300);
                if (orf_len >= MIN_ORF_BP) {
                    candidates.push_back({node, dist + offset, offset, orf_len});
                }
            }
        }
        
        if (dist < max_dist) {
            uint64_t edges[4];
            int edge_count;
            
            // Direction-aware edge selection
            if (direction == SearchDirection::DOWNSTREAM) {
                edge_count = sdbg_.OutgoingEdges(node, edges);
            } else {
                edge_count = sdbg_.IncomingEdges(node, edges);
            }
            
            for (int i = 0; i < edge_count; ++i) {
                if (sdbg_.IsValidEdge(edges[i]) && !visited.count(edges[i])
                    && (cycle_nodes_.empty() || !cycle_nodes_[edges[i]])) {
                    q.push({edges[i], dist + 1});
                    visited.insert(edges[i]);
                }
            }
        }
    }
    return candidates;
}

// PROGRESSIVE DEPTH: Quick shallow scan to reject non-promising candidates
// Returns best normalized score found in the first shallow window.
// Window = max(120 bp, 60% of profile length in bp) so that proteins whose
// conserved domain sits past position ~40 aa are not prematurely pruned.
double CasWorkflow::ScoreStartNodeShallow(
    uint64_t start_node, int start_codon_offset,
    size_t profile_index, SearchDirection direction) {
    const Profile* hmm = GetProfile(profile_index);
    if (!hmm) return 0.0;

    // Scale with profile length: scan at least 60% of the model before deciding
    const int shallow_depth = std::max(120, hmm->GetLength() * 3 * 3 / 5);  // aa→bp, 60%

    CasGeneDetector detector(sdbg_, hmm);
    auto paths = detector.BeamSearchAminoAcids(
        start_node, params_.BEAM_WIDTH, shallow_depth, start_codon_offset,
        direction);
    
    double best_score = 0.0;
    for (const auto& path : paths) {
        if (path.amino_acids.empty()) continue;
        // Use partial normalized score (score / positions_matched)
        // Even incomplete paths give a signal
        if (path.normalized_score > best_score) {
            best_score = path.normalized_score;
        }
    }
    return best_score;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithProfile(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    size_t profile_index, SearchDirection direction,
    bool skip_distance_check) {
    DetectedCasGene result;
    result.start_node = start_node;
    result.distance_from_repeat = distance_from_repeat;
    result.profile_name = profiles_[profile_index].filename;
    result.profile_index = profile_index;
    
    // CRITICAL: Check if profile is compatible with available space (UPSTREAM filtering)
    // Skip when caller has already validated (skip_distance_check=true)
    if (!skip_distance_check && !IsProfileCompatibleWithDistance(profiles_[profile_index], distance_from_repeat, direction)) {
        return result;  // Profile too long for available space
    }
    
    // NOTE: ORF length check removed - caller should pre-compute and filter
    // This avoids expensive BFS inside the hot parallel loop
    
    // Use pre-loaded profile (thread-safe read-only access)
    const Profile* hmm = GetProfile(profile_index);
    if (!hmm) {
        return result;  // Profile failed to load
    }
    
    CasGeneDetector detector(sdbg_, hmm);
    auto paths = detector.BeamSearchAminoAcids(
        start_node, params_.BEAM_WIDTH, profiles_[profile_index].max_bp, start_codon_offset,
        direction);
    
    for (const auto& path : paths) {
        if (path.amino_acids.empty() || !path.is_complete) continue;
        if (static_cast<int>(path.dna_sequence.size()) > profiles_[profile_index].max_bp) continue;
        if (static_cast<int>(path.dna_sequence.size()) < profiles_[profile_index].min_bp) continue;

        if (path.normalized_score > result.normalized_score) {
            result.bit_score = path.total_score;
            result.normalized_score = path.normalized_score;
            result.gene_length = static_cast<int>(path.dna_sequence.size());
            result.end_node = path.node_path.back();
            result.node_path = path.node_path;
            result.is_complete = true;
            // Efficient amino acid string building with pre-allocation
            result.amino_acids.clear();
            result.amino_acids.reserve(path.amino_acids.size());
            for (const auto& aa : path.amino_acids) result.amino_acids += aa;
        }
    }
    
    return result;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithAllProfiles(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    SearchDirection direction) {

    DetectedCasGene best;

    // Test ALL profiles - no shortcuts
    std::vector<DetectedCasGene> results(profiles_.size());
    for (size_t i = 0; i < profiles_.size(); ++i) {
        results[i] = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, i, direction);
        results[i].gene_family = ExtractGeneFamily(profiles_[i].filename);
    }

    // Deduplicate: keep best score per family
    std::map<std::string, DetectedCasGene> best_per_family;
    for (const auto& r : results) {
        if (r.is_complete) {
            auto it = best_per_family.find(r.gene_family);
            if (it == best_per_family.end() || r.normalized_score > it->second.normalized_score) {
                best_per_family[r.gene_family] = r;
            }
        }
    }

    for (const auto& [family, gene] : best_per_family) {
        if (gene.normalized_score > best.normalized_score) {
            best = gene;
        }
    }
    return best;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithExpectedFamilies(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    const std::set<std::string>& expected_families, SearchDirection direction) {

    DetectedCasGene best;
    
    // Only test profiles matching expected families
    std::vector<std::pair<size_t, std::string>> to_test;
    for (const auto& family : expected_families) {
        auto it = family_to_profiles_.find(family);
        if (it != family_to_profiles_.end()) {
            for (size_t idx : it->second) {
                to_test.push_back({idx, family});
            }
        }
    }
    
    if (to_test.empty()) return best;
    
    std::vector<DetectedCasGene> results(to_test.size());
    for (size_t i = 0; i < to_test.size(); ++i) {
        results[i] = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, 
            to_test[i].first, direction);
        results[i].gene_family = to_test[i].second;
    }

    // Deduplicate: keep best score per family
    std::map<std::string, DetectedCasGene> best_per_family;
    for (const auto& r : results) {
        if (r.is_complete) {
            auto it = best_per_family.find(r.gene_family);
            if (it == best_per_family.end() || r.normalized_score > it->second.normalized_score) {
                best_per_family[r.gene_family] = r;
            }
        }
    }

    for (const auto& [family, gene] : best_per_family) {
        if (gene.normalized_score > best.normalized_score) {
            best = gene;
        }
    }
    return best;
}

DetectedCasGene CasWorkflow::ScoreStartNodeExploratory(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    const std::set<std::string>& detected_families, SearchDirection direction) {

    DetectedCasGene best;
    
    // Test all profiles EXCEPT already-detected families
    std::vector<std::pair<size_t, std::string>> to_test;
    for (const auto& [family, indices] : family_to_profiles_) {
        if (!detected_families.count(family)) {
            for (size_t idx : indices) {
                to_test.push_back({idx, family});
            }
        }
    }
    
    if (to_test.empty()) return best;
    
    std::vector<DetectedCasGene> results(to_test.size());
    for (size_t i = 0; i < to_test.size(); ++i) {
        results[i] = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, 
            to_test[i].first, direction);
        results[i].gene_family = to_test[i].second;
        results[i].is_putative = true;  // Exploratory finds are putative
    }

    // Deduplicate: keep best score per family
    std::map<std::string, DetectedCasGene> best_per_family;
    for (const auto& r : results) {
        if (r.is_complete) {
            auto it = best_per_family.find(r.gene_family);
            if (it == best_per_family.end() || r.normalized_score > it->second.normalized_score) {
                best_per_family[r.gene_family] = r;
            }
        }
    }

    for (const auto& [family, gene] : best_per_family) {
        if (gene.normalized_score > best.normalized_score) {
            best = gene;
        }
    }
    return best;
}

// =============================================================================
// MULTI-PROFILE SCORING: Score ONE candidate against ALL specified profiles
// in a single graph traversal. 100x fewer graph traversals!
// =============================================================================
std::vector<std::pair<size_t, DetectedCasGene>> CasWorkflow::ScoreStartNodeMultiProfile(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    const std::vector<std::pair<size_t, std::string>>& profile_family_pairs,
    SearchDirection direction) {
    
    std::vector<std::pair<size_t, DetectedCasGene>> results;
    if (profile_family_pairs.empty()) return results;
    
    // Build profile list for BeamSearchMultiProfile
    std::vector<std::pair<const Profile*, size_t>> profiles_for_beam;
    profiles_for_beam.reserve(profile_family_pairs.size());
    
    for (const auto& [profile_idx, family] : profile_family_pairs) {
        const Profile* p = GetProfile(profile_idx);
        if (p) {
            profiles_for_beam.push_back({p, profile_idx});
        }
    }
    
    if (profiles_for_beam.empty()) return results;
    
    // Find max profile length for depth calculation
    int max_len = 0;
    for (const auto& [p, idx] : profiles_for_beam) {
        max_len = std::max(max_len, p->GetLength());
    }
    int max_depth = static_cast<int>(max_len * 1.3 * 3);  // AA to bp, with margin
    
    // Create detector and run multi-profile beam search
    CasGeneDetector detector(sdbg_, nullptr);  // No single profile
    auto multi_results = detector.BeamSearchMultiProfile(
        start_node, params_.BEAM_WIDTH, max_depth, start_codon_offset,
        profiles_for_beam, params_.MIN_NORMALIZED_SCORE * 0.5);  // Lower threshold for beam
    
    // Convert MultiProfileResult to DetectedCasGene for each profile
    for (const auto& mr : multi_results) {
        for (size_t i = 0; i < profile_family_pairs.size(); ++i) {
            size_t profile_idx = profile_family_pairs[i].first;
            const std::string& family = profile_family_pairs[i].second;
            
            if (i >= mr.profile_scores.size()) continue;
            
            auto [total_score, normalized_score, hmm_pos] = mr.profile_scores[i];
            
            if (normalized_score >= params_.MIN_NORMALIZED_SCORE) {
                DetectedCasGene gene;
                gene.profile_name = profiles_[profile_idx].filename;
                gene.profile_index = profile_idx;
                gene.gene_family = family;
                gene.start_node = start_node;
                gene.distance_from_repeat = distance_from_repeat;
                gene.gene_length = mr.amino_acids.size() * 3;
                gene.normalized_score = normalized_score;
                gene.bit_score = total_score;
                gene.amino_acids = mr.amino_acids;
                gene.node_path = mr.node_path;
                gene.is_complete = mr.is_complete;
                gene.is_putative = false;
                
                results.push_back({profile_idx, gene});
            }
        }
    }
    
    return results;
}

DetectedCasGene CasWorkflow::ScoreStartNodeWithFirstGeneProfiles(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    SearchDirection direction) {
    
    DetectedCasGene best;
    
    // CRITICAL: If no first-gene profiles, detection fails silently
    if (first_gene_profile_indices_.empty()) {
        return best;  // Empty result - no profiles to test!
    }
    
    // Test all first-gene profiles (they're already filtered to ~22 profiles)
    for (size_t profile_idx : first_gene_profile_indices_) {
        auto result = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, profile_idx, direction);
        
        if (result.is_complete) {
            result.gene_family = ExtractGeneFamily(profiles_[profile_idx].filename);
            
            if (result.normalized_score > best.normalized_score) {
                best = result;
            }
        }
    }
    
    return best;
}

// SEPARATE FUNCTION FOR DOWNSTREAM - uses shallow pre-filtering to handle many candidates
DetectedCasGene CasWorkflow::ScoreStartNodeWithFirstGeneProfilesDownstream(
    uint64_t start_node, int distance_from_repeat, int start_codon_offset,
    SearchDirection direction) {
    
    DetectedCasGene best;
    
    if (first_gene_profile_indices_.empty()) {
        return best;
    }
    
    const double shallow_threshold = tuning_.shallow_threshold;
    
    // For DOWNSTREAM: pre-filter with shallow scan before full HMM
    for (size_t profile_idx : first_gene_profile_indices_) {
        // Quick shallow scan first (skip for short profiles where N-terminus may lack signal)
        if (profiles_[profile_idx].leng >= 200) {
            double shallow_score = ScoreStartNodeShallow(
                start_node, start_codon_offset, profile_idx, direction);
            if (shallow_score < shallow_threshold) {
                continue;  // Reject - no signal
            }
        }
        
        // Passed shallow filter - do full HMM scan
        auto result = ScoreStartNodeWithProfile(
            start_node, distance_from_repeat, start_codon_offset, profile_idx, direction);
        
        if (result.is_complete) {
            result.gene_family = ExtractGeneFamily(profiles_[profile_idx].filename);
            
            if (result.normalized_score > best.normalized_score) {
                best = result;
            }
        }
    }
    
    return best;
}

std::vector<size_t> CasWorkflow::GetProfileIndicesForType(const std::string& type_class) const {
    std::vector<size_t> indices;
    std::set<size_t> seen;
    
    // Find the rule for this type
    for (const auto& rule : type_rules_) {
        if (rule.type_class == type_class) {
            // Helper to resolve a rule family name to all matching profile indices.
            // Rule entries are families (e.g. cas5), not exact profile stems (e.g. cas5_I-E_2).
            auto resolveGene = [&](const std::string& gene_entry) {
                if (gene_entry.find('|') != std::string::npos) {
                    std::stringstream ss(gene_entry);
                    std::string alt;
                    while (std::getline(ss, alt, '|')) {
                        std::string family = ExtractGeneFamily(alt);
                        auto it = family_to_profiles_.find(family);
                        if (it != family_to_profiles_.end()) {
                            for (size_t idx : it->second) {
                                if (seen.insert(idx).second) {
                                    indices.push_back(idx);
                                }
                            }
                        }
                    }
                } else {
                    std::string family = ExtractGeneFamily(gene_entry);
                    auto it = family_to_profiles_.find(family);
                    if (it != family_to_profiles_.end()) {
                        for (size_t idx : it->second) {
                            if (seen.insert(idx).second) {
                                indices.push_back(idx);
                            }
                        }
                    }
                }
            };
            
            // Get all profile indices for genes in this rule
            for (const auto& gene_name : rule.mandatory_genes) {
                resolveGene(gene_name);
            }
            for (const auto& gene_name : rule.accessory_genes) {
                resolveGene(gene_name);
            }
            break;
        }
    }
    
    return indices;
}

// Score how well a set of detected gene families matches a candidate type
// Returns: pair<interference_matches, total_matches>
std::pair<int, int> CasWorkflow::ScoreTypeMatch(
    size_t type_idx, const std::set<std::string>& detected_families) const {
    
    if (type_idx >= type_rules_.size()) return {0, 0};
    const auto& rule = type_rules_[type_idx];
    
    int interf_matches = rule.CountMatchedInterference(detected_families);
    int total_matches = rule.CountMatchedTotal(detected_families);
    
    return {interf_matches, total_matches};
}

// Find all candidate types that could match a set of detected gene families
std::vector<size_t> CasWorkflow::FindCandidateTypes(
    const std::set<std::string>& detected_families) const {
    
    std::vector<size_t> candidates;
    
    for (size_t i = 0; i < type_rules_.size(); ++i) {
        const auto& rule = type_rules_[i];
        int interf_matches = rule.CountMatchedInterference(detected_families);
        
        // Must have at least one interference gene match to be a candidate
        if (interf_matches > 0) {
            candidates.push_back(i);
        }
    }
    
    return candidates;
}

// Narrow candidate types based on a newly detected gene family.
// Multi-model families (cas1/2/3/4/6, regulatory) appear in almost every type,
// so detecting them provides no discriminatory power — skip narrowing for them.
std::vector<size_t> CasWorkflow::NarrowCandidateTypesByFamily(
    const std::vector<size_t>& candidate_type_indices,
    const std::string& detected_family) const {

    // Do not narrow on ubiquitous multi-model genes
    if (MULTI_MODEL_FAMILIES.count(detected_family)) {
        return candidate_type_indices;
    }

    std::vector<size_t> remaining;
    for (size_t type_idx : candidate_type_indices) {
        if (type_idx >= type_rules_.size()) continue;
        const auto& rule = type_rules_[type_idx];
        if (rule.MatchesInterferenceFamily(detected_family) ||
            rule.MatchesAdaptationFamily(detected_family) ||
            rule.MatchesAccessoryFamily(detected_family)) {
            remaining.push_back(type_idx);
        }
    }
    return remaining;
}

// Legacy method - detects only downstream for backwards compatibility
std::vector<DetectedCasGene> CasWorkflow::DetectCasGenes(uint64_t repeat_node) {
    std::vector<DetectedCasGene> genes;

    auto candidates = FindStartCodonCandidates(
        repeat_node, params_.FIRST_GENE_MIN_DIST, params_.FIRST_GENE_MAX_DIST,
        params_.MAX_START_CANDIDATES, SearchDirection::DOWNSTREAM);

    if (candidates.empty()) return genes;

    std::vector<DetectedCasGene> results(candidates.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < candidates.size(); ++i) {
        const auto& c = candidates[i];
        results[i] = ScoreStartNodeWithAllProfiles(c.node, c.distance, c.offset, 
                                                   SearchDirection::DOWNSTREAM);
    }

    // Deduplicate: keep best per family
    std::map<std::string, DetectedCasGene> best_per_family;
    for (const auto& r : results) {
        if (r.normalized_score >= (r.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(r.profile_index))) {
            auto it = best_per_family.find(r.gene_family);
            if (it == best_per_family.end() || r.normalized_score > it->second.normalized_score) {
                best_per_family[r.gene_family] = r;
            }
        }
    }

    // Return best overall
    DetectedCasGene best;
    for (const auto& [family, gene] : best_per_family) {
        if (gene.normalized_score > best.normalized_score) {
            best = gene;
        }
    }

    if (best.is_complete) {
        genes.push_back(best);
    }
    return genes;
}

CasCassette CasWorkflow::DetectCassette(uint64_t repeat_node, SearchDirection direction) {
    std::cout << __func__ << "\n";
    CasCassette cassette;
    cassette.repeat_node = repeat_node;
    cassette.direction = direction;
    cassette.stop_reason_code.clear();

    const int max_cassette_bp = tuning_.max_cassette_bp;
    const int intergenic_min = tuning_.intergenic_min_bp;   // Allow overlapping genes
    const int intergenic_max = tuning_.intergenic_max_bp;   // Max intergenic gap
    const int max_genes = tuning_.max_genes;
    const int profile_window_count = tuning_.profile_window_count;

    std::set<std::string> detected_families;  // For deduplication
    
    DetectedCasGene gene1;
    size_t original_profile_index = SIZE_MAX;
    constexpr size_t FIRST_GENE_TOP_K = 10;
    
    // STEP 1: First gene search (profile-parallel)
    // IMPORTANT: avoid a single huge serial BFS over [min..5000], which bottlenecks on one core.
    // Each profile searches only its own expected distance window in parallel.
    int global_search_min = params_.FIRST_GENE_MIN_DIST;
    int global_search_max = tuning_.first_gene_search_max_bp;

    std::cout << " (" << first_gene_profile_indices_.size() << " first-gene profiles)..." << std::flush;

    // =========================================================================
    // OPTIMIZATION: Single BFS to find ALL candidates, then filter per profile
    // =========================================================================
    auto all_first_gene_candidates = FindStartCodonCandidates(
        repeat_node, global_search_min, global_search_max, tuning_.bfs_max_candidates, direction);
    
    // Sort by distance for efficient filtering
    std::sort(all_first_gene_candidates.begin(), all_first_gene_candidates.end(),
        [](const StartCodonCandidate& a, const StartCodonCandidate& b) {
            return a.distance < b.distance;
        });

    std::vector<DetectedCasGene> best_by_profile(first_gene_profile_indices_.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t pi = 0; pi < first_gene_profile_indices_.size(); ++pi) {
        size_t profile_idx = first_gene_profile_indices_[pi];
        const auto& profile = profiles_[profile_idx];

        // FIRST GENE: Distance from repeat is INDEPENDENT of gene length!
        // Don't filter by profile length - the gene can be anywhere in the search range.
        // Just use the global range and let HMM scoring decide.
        int allowed_min = global_search_min;
        int allowed_max = global_search_max;

        // Use more bins for first-gene since we're searching a wide range (0-5000bp)
        // 20 bins = ~250bp per bin, giving reasonable coverage
        const int first_gene_bins = tuning_.first_gene_bins;

        // Filter from pre-computed candidates (no BFS here!)
        auto profile_candidates = BuildBinnedCandidateWindow(
            all_first_gene_candidates, allowed_min, allowed_max, first_gene_bins);

        DetectedCasGene local_best;
        for (const auto* cptr : profile_candidates) {
            const auto& c = *cptr;

            auto result = ScoreStartNodeWithProfile(
                c.node, c.distance, c.offset, profile_idx, direction, true);

            if (result.normalized_score >= (result.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(result.profile_index))) {
                result.gene_family = ExtractGeneFamily(profile.filename);
                if (result.normalized_score > local_best.normalized_score) {
                    local_best = result;
                }
            }
        }
        best_by_profile[pi] = local_best;
    }

    std::vector<DetectedCasGene> first_gene_hypotheses;
    first_gene_hypotheses.reserve(best_by_profile.size());
    for (const auto& hit : best_by_profile) {
        if (hit.normalized_score >= (hit.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(hit.profile_index))) {
            first_gene_hypotheses.push_back(hit);
        }
    }

    std::sort(first_gene_hypotheses.begin(), first_gene_hypotheses.end(),
        [](const DetectedCasGene& a, const DetectedCasGene& b) {
            return a.normalized_score > b.normalized_score;
        });

    if (first_gene_hypotheses.size() > FIRST_GENE_TOP_K) {
        first_gene_hypotheses.resize(FIRST_GENE_TOP_K);
    }

    std::cout << "  DEBUG: first-gene hypotheses after scoring="
              << first_gene_hypotheses.size() << std::endl;

    // Resolve TOP-K hypotheses after family refinement; keep source profile index.
    for (auto hypothesis : first_gene_hypotheses) {
        size_t hypothesis_source_profile_index = hypothesis.profile_index;

        auto family_it = family_to_profiles_.find(hypothesis.gene_family);
        if (family_it != family_to_profiles_.end() && family_it->second.size() > 1) {
            int start_offset = GetStartCodonOffset(hypothesis.start_node);
            if (start_offset < 0) start_offset = 0;

            std::vector<DetectedCasGene> family_results(family_it->second.size());

            #pragma omp parallel for schedule(dynamic)
            for (size_t f = 0; f < family_it->second.size(); ++f) {
                size_t profile_idx = family_it->second[f];
                auto fresult = ScoreStartNodeWithProfile(
                    hypothesis.start_node, hypothesis.distance_from_repeat, start_offset,
                    profile_idx, direction, true);
                if (fresult.is_complete || fresult.normalized_score >= PartialScoreThreshold(fresult.profile_index)) {
                    fresult.gene_family = hypothesis.gene_family;
                    family_results[f] = fresult;
                }
            }

            for (const auto& fresult : family_results) {
                if ((fresult.is_complete || fresult.normalized_score >= PartialScoreThreshold(fresult.profile_index)) && fresult.normalized_score > hypothesis.normalized_score) {
                    hypothesis = fresult;
                }
            }
        }

        if (!gene1.is_complete || hypothesis.normalized_score > gene1.normalized_score) {
            gene1 = hypothesis;
            original_profile_index = hypothesis_source_profile_index;
        }
    }

    if (gene1.normalized_score < (gene1.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(gene1.profile_index))) {
        cassette.stop_reason_code = "NO_NEXT_GENE";
        return cassette;
    }

    // Post-output length validation: verify amino acid length vs profile bounds
    if (gene1.profile_index < profiles_.size()) {
        int aa_len = static_cast<int>(gene1.amino_acids.size());
        const auto& prof = profiles_[gene1.profile_index];
        if (aa_len < prof.min_aa || aa_len > prof.max_aa) {
            std::cout << "  WARNING: Gene1 " << gene1.gene_family
                      << " AA length " << aa_len << " outside profile bounds ["
                      << prof.min_aa << "," << prof.max_aa << "] - rejected" << std::endl;
            cassette.stop_reason_code = "NO_NEXT_GENE";
            return cassette;
        }
    }

    // =========================================================================
    // Get ALL candidate type indices that have this first gene
    // Use original_profile_index since family refinement may have changed it
    // =========================================================================
    std::vector<size_t> candidate_type_indices;
    
    auto type_it = first_gene_to_types_.find(original_profile_index);
    std::cout << "  DEBUG: original_profile_index=" << original_profile_index 
              << " profile=" << profiles_[original_profile_index].filename << std::endl;
    if (type_it != first_gene_to_types_.end()) {
        std::cout << "  DEBUG: Found " << type_it->second.size() << " types for this profile: ";
        for (const std::string& type_name : type_it->second) {
            std::cout << type_name << " ";
            for (size_t i = 0; i < type_rules_.size(); ++i) {
                if (type_rules_[i].type_class == type_name) {
                    candidate_type_indices.push_back(i);
                    break;
                }
            }
        }
        std::cout << std::endl;
    } else {
        std::cout << "  DEBUG: NO types found in first_gene_to_types_ for profile index " 
                  << original_profile_index << std::endl;
    }

    cassette.genes.push_back(gene1);
    cassette.total_distance_bp = gene1.distance_from_repeat + gene1.gene_length;
    detected_families.insert(gene1.gene_family);
    

    // =========================================================================
    // SUBSEQUENT GENES: Position-aware search
    // For UPSTREAM: only check genes BEFORE current position in rule
    // For DOWNSTREAM: only check genes AFTER current position in rule
    // =========================================================================
    for (int iter = 0; iter < max_genes && cassette.total_distance_bp < max_cassette_bp; ++iter) {
        const DetectedCasGene& prev_gene = cassette.genes.back();
        
        // =====================================================================
        // PROFILE-LENGTH-TARGETED SEARCH (PLTS) - Profile-Major
        // =====================================================================
        // For each profile:
        //   1. Calculate WHERE its START codon should be: D = profile_length + gap
        //   2. Find STARTs only in [D - gap_max, D + gap_tolerance] (~300bp window)
        //   3. Shallow scan - if promising, do full scan
        //   4. If score > threshold, extend ±25% to find optimal
        //
        // UPSTREAM: D = next_gene_length + gap, gap ∈ [-24, 300]
        //   So for profile of length L: search [L - 24, L + 300]
        //
        // DOWNSTREAM: D = prev_gene_length + gap
        //   Search window is just [prev_len - 24, prev_len + 300]
        // =====================================================================
        
        uint64_t search_origin = prev_gene.start_node;
        int prev_gene_len = prev_gene.gene_length;
        
        
        DetectedCasGene next_gene;
        
        // ---------------------------------------------------------------------
        // Collect expected families from candidate types
        // ---------------------------------------------------------------------
        std::set<std::string> expected_families;
        if (!candidate_type_indices.empty()) {
            for (size_t type_idx : candidate_type_indices) {
                const auto& rule = type_rules_[type_idx];
                for (const auto& alt_set : rule.mandatory_families) {
                    for (const auto& fam : alt_set) {
                        if (detected_families.count(fam) == 0) {
                            expected_families.insert(fam);
                        }
                    }
                }
                for (const auto& alt_set : rule.accessory_families) {
                    for (const auto& fam : alt_set) {
                        if (detected_families.count(fam) == 0) {
                            expected_families.insert(fam);
                        }
                    }
                }
            }
        }
        
        
        // If no expected families from type rules, skip to Phase 2 exploratory
        bool skip_phase1 = expected_families.empty();
        
        std::cout << "  DEBUG iter " << iter << ": prev_gene=" << prev_gene.gene_family 
                  << " dist=" << prev_gene.distance_from_repeat << " len=" << prev_gene.gene_length
                  << " | candidate_types=" << candidate_type_indices.size()
                  << " expected_families=" << expected_families.size()
                  << " skip_phase1=" << skip_phase1 << std::endl;
        
        // ---------------------------------------------------------------------
        // PROFILE-MAJOR SEARCH (profile-parallel candidate discovery)
        // ---------------------------------------------------------------------
        std::vector<std::pair<size_t, std::string>> representative_profiles;
        
        for (const auto& fam : expected_families) {
            auto it = family_to_profiles_.find(fam);
            if (it != family_to_profiles_.end() && !it->second.empty()) {
                // Include ALL profiles for the family, not just the first
                for (size_t idx : it->second) {
                    representative_profiles.push_back({idx, fam});
                }
            }
        }

        // =========================================================================
        // PHASE 1: Type-guided search (skip if no expected families)
        // GEOMETRY:
        //   UPSTREAM: dist(prev_start -> next_start) = next_gene_len + gap
        //   DOWNSTREAM: dist(prev_start -> next_start) = prev_gene_len + gap
        // =========================================================================
        
        // Declare outside of conditional block so they're accessible later
        DetectedCasGene best_gene;
        
        if (!skip_phase1 && !representative_profiles.empty()) {
            int phase1_search_min, phase1_search_max;
            
            if (direction == SearchDirection::DOWNSTREAM) {
                // DOWNSTREAM: search range is constant (depends only on prev_gene_len)
                phase1_search_min = prev_gene_len + intergenic_min;
                phase1_search_max = prev_gene_len + intergenic_max;
            } else {
                // UPSTREAM: search range depends on profile lengths
                // dist = next_gene_len + gap
                //   next_gene_len ∈ [min_bp, max_bp] (already -15%/+25% of LENG*3)
                //   gap ∈ [INTERGENIC_MIN, INTERGENIC_MAX]
                phase1_search_min = std::numeric_limits<int>::max();
                phase1_search_max = 0;
                
                for (const auto& [profile_idx, family] : representative_profiles) {
                    const auto& profile = profiles_[profile_idx];
                    int prof_min = profile.min_bp + intergenic_min;
                    int prof_max = profile.max_bp + intergenic_max;
                    phase1_search_min = std::min(phase1_search_min, prof_min);
                    phase1_search_max = std::max(phase1_search_max, prof_max);
                }
            }
        
        // Single BFS to find ALL candidates in the widest range
        auto all_phase1_candidates = FindStartCodonCandidatesFromNode(
            search_origin, phase1_search_min, phase1_search_max, tuning_.bfs_max_candidates, direction);
        
        std::cout << "  DEBUG Phase1: profiles=" << representative_profiles.size()
                  << " search_range=[" << phase1_search_min << "," << phase1_search_max << "]"
                  << " BFS found " << all_phase1_candidates.size() << " candidates" << std::endl;
        
            // Sort by distance for efficient filtering
            std::sort(all_phase1_candidates.begin(), all_phase1_candidates.end(),
                [](const StartCodonCandidate& a, const StartCodonCandidate& b) {
                    return a.distance < b.distance;
                });

            // For each profile, filter to matching distance range and score
            std::atomic<bool> found_good_gene{false};
            std::mutex best_mutex;
            
            #pragma omp parallel for schedule(dynamic)
            for (size_t pi = 0; pi < representative_profiles.size(); ++pi) {
            if (found_good_gene.load(std::memory_order_relaxed)) continue;
            
            const auto& [profile_idx, family] = representative_profiles[pi];
            const auto& profile = profiles_[profile_idx];
            
            // Direction-aware distance range for THIS profile
            // CRITICAL GEOMETRY:
            //   UPSTREAM: dist = next_gene_len + gap
            //     - next_gene_len ∈ [min_bp, max_bp] (already -15%/+25% of LENG*3)
            //     - gap ∈ [INTERGENIC_MIN, INTERGENIC_MAX]
            //     - So dist ∈ [min_bp + INTERGENIC_MIN, max_bp + INTERGENIC_MAX]
            //   DOWNSTREAM: dist = prev_gene_len + gap → independent of profile_bp!
            int allowed_min, allowed_max;
            if (direction == SearchDirection::UPSTREAM) {
                // Use min_bp/max_bp directly (already encode -15%/+25% tolerance)
                allowed_min = profile.min_bp + intergenic_min;
                // Maximum: longest tolerated gene + maximum gap
                allowed_max = profile.max_bp + intergenic_max;
            } else {
                // DOWNSTREAM: distance is fixed regardless of which profile we're testing
                allowed_min = prev_gene_len + intergenic_min;
                allowed_max = prev_gene_len + intergenic_max;
            }

            // Filter from pre-computed candidates (no BFS here!)
            auto profile_candidates = BuildBinnedCandidateWindow(
                all_phase1_candidates, allowed_min, allowed_max, profile_window_count);
            
            // Test ALL candidates that pass shallow threshold, not just the best one
            // This avoids missing genes when a wrong START has higher shallow score
            for (const auto* cptr : profile_candidates) {
                if (found_good_gene.load(std::memory_order_relaxed)) break;
                
                bool skip_shallow = (profiles_[profile_idx].leng < 200);
                double shallow_score = skip_shallow ? tuning_.shallow_threshold :
                    ScoreStartNodeShallow(cptr->node, cptr->offset, profile_idx, direction);
                
                if (shallow_score >= tuning_.shallow_threshold) {
                    const auto& c = *cptr;
                    int global_dist = prev_gene.distance_from_repeat + c.distance;
                    
                    auto result = ScoreStartNodeWithProfile(
                        c.node, global_dist, c.offset, profile_idx, direction, true);
                    
                    if (result.normalized_score >= (result.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(result.profile_index))) {
                        result.gene_family = family;
                        std::lock_guard<std::mutex> lock(best_mutex);
                        if (result.normalized_score > best_gene.normalized_score) {
                            best_gene = result;
                            // Only early exit on very high confidence hits
                            if (result.normalized_score >= tuning_.phase1_high_confidence) {
                                found_good_gene.store(true, std::memory_order_relaxed);
                            }
                        }
                    }
                }
            }
        }  // end omp parallel for
        }  // end if (!skip_phase1)
        
        std::cout << "  DEBUG Phase1 result: best_gene.is_complete=" << best_gene.is_complete
                  << " score=" << best_gene.normalized_score << std::endl;
        
        if (best_gene.is_complete || best_gene.normalized_score >= PartialScoreThreshold(best_gene.profile_index)) {
            next_gene = best_gene;
        }
        
        // FAMILY REFINEMENT if we found something
        if (next_gene.normalized_score >= (next_gene.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(next_gene.profile_index))) {
            
            auto family_it = family_to_profiles_.find(next_gene.gene_family);
            if (family_it != family_to_profiles_.end() && family_it->second.size() > 1) {
                
                // Get the actual start codon offset within the node
                int start_offset = GetStartCodonOffset(next_gene.start_node);
                if (start_offset < 0) start_offset = 0;  // Fallback
                
                // Use the already-found start node directly (no re-search needed)
                std::vector<DetectedCasGene> family_results(family_it->second.size());
                
                #pragma omp parallel for schedule(dynamic)
                for (size_t f = 0; f < family_it->second.size(); ++f) {
                    size_t pidx = family_it->second[f];
                    // Score at the SAME start node with correct offset
                    auto result = ScoreStartNodeWithProfile(
                        next_gene.start_node, next_gene.distance_from_repeat, 
                        start_offset,
                        pidx, direction, true);
                    if (result.is_complete || result.normalized_score >= PartialScoreThreshold(result.profile_index)) {
                        result.gene_family = next_gene.gene_family;
                        family_results[f] = result;
                    }
                }
                
                for (const auto& result : family_results) {
                    if ((result.is_complete || result.normalized_score >= PartialScoreThreshold(result.profile_index)) && result.normalized_score > next_gene.normalized_score) {
                        next_gene = result;
                    }
                }
            }
        }

        // ---------------------------------------------------------------------
        // PHASE 2: Exploratory search (if Phase 1 fails OR was skipped)
        // Uses profile-major approach - tests ALL remaining families
        // ---------------------------------------------------------------------
        if (skip_phase1 || (params_.allow_exploratory &&
            (next_gene.normalized_score < (next_gene.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(next_gene.profile_index))))) {
            
            std::cout << "  DEBUG: Entering Phase 2 (skip_phase1=" << skip_phase1 
                      << " allow_exploratory=" << params_.allow_exploratory << ")" << std::endl;
            
            // Collect profiles not yet detected
            std::vector<std::pair<size_t, std::string>> profiles_to_explore;
            for (const auto& [family, indices] : family_to_profiles_) {
                if (detected_families.count(family)) continue;  // Already detected
                for (size_t idx : indices) {
                    profiles_to_explore.push_back({idx, family});
                }
            }
            
            // =========================================================================
            // OPTIMIZATION: Compute global search range and find candidates ONCE
            // =========================================================================
            int global_search_min, global_search_max;
            
            if (direction == SearchDirection::DOWNSTREAM) {
                // DOWNSTREAM: search range is constant (depends only on prev_gene_len)
                global_search_min = prev_gene_len + intergenic_min;
                global_search_max = prev_gene_len + intergenic_max;
            } else {
                // UPSTREAM: search range depends on profile lengths
                // dist = next_gene_len + gap
                //   next_gene_len ∈ [min_bp, max_bp] (already -15%/+25% of LENG*3)
                //   gap ∈ [INTERGENIC_MIN, INTERGENIC_MAX]
                global_search_min = std::numeric_limits<int>::max();
                global_search_max = 0;
                
                for (const auto& [profile_idx, family] : profiles_to_explore) {
                    const auto& profile = profiles_[profile_idx];
                    int search_min = profile.min_bp + intergenic_min;
                    int search_max = profile.max_bp + intergenic_max;
                    global_search_min = std::min(global_search_min, search_min);
                    global_search_max = std::max(global_search_max, search_max);
                }
            }
            
            // Single BFS to find ALL candidates in the widest range
            auto all_candidates = FindStartCodonCandidatesFromNode(
                search_origin, global_search_min, global_search_max, 
                tuning_.bfs_max_candidates, direction);
            
            std::cout << "  DEBUG Phase2: profiles=" << profiles_to_explore.size()
                      << " search_range=[" << global_search_min << "," << global_search_max << "]"
                      << " BFS found " << all_candidates.size() << " candidates" << std::endl;
            
            // Sort by distance for BuildBinnedCandidateWindow's binary search
            std::sort(all_candidates.begin(), all_candidates.end(),
                [](const StartCodonCandidate& a, const StartCodonCandidate& b) {
                    return a.distance < b.distance;
                });
            
            // Profile-major search: for each profile, filter candidates by distance
            DetectedCasGene best_p2;
            std::atomic<bool> found_high_confidence_p2{false};
            std::mutex p2_mutex;
            
            #pragma omp parallel for schedule(dynamic)
            for (int p = 0; p < static_cast<int>(profiles_to_explore.size()); ++p) {
                if (found_high_confidence_p2.load(std::memory_order_relaxed)) continue;
                const auto& [profile_idx, family] = profiles_to_explore[p];
                const auto& profile = profiles_[profile_idx];
                
                // Calculate allowed distance range for this profile
                // CRITICAL GEOMETRY:
                //   UPSTREAM: dist = next_gene_len + gap
                //     - next_gene_len ∈ [min_bp, max_bp] (already -15%/+25% of LENG*3)
                //     - gap ∈ [INTERGENIC_MIN, INTERGENIC_MAX]
                //   DOWNSTREAM: dist = prev_gene_len + gap → independent of profile!
                int allowed_min, allowed_max;
                if (direction == SearchDirection::UPSTREAM) {
                    allowed_min = profile.min_bp + intergenic_min;
                    allowed_max = profile.max_bp + intergenic_max;
                } else {
                    // DOWNSTREAM: distance is fixed regardless of which profile we're testing
                    allowed_min = prev_gene_len + intergenic_min;
                    allowed_max = prev_gene_len + intergenic_max;
                }
                
                // Filter from pre-computed candidates (no BFS here!)
                auto profile_candidates = BuildBinnedCandidateWindow(
                    all_candidates, allowed_min, allowed_max, profile_window_count);
                if (profile_candidates.empty()) continue;
                
                // Test this profile against candidates
                DetectedCasGene local_best_p2;
                for (const auto* cptr : profile_candidates) {
                    if (found_high_confidence_p2.load(std::memory_order_relaxed)) break;

                    const auto& c = *cptr;

                    // Skip ORF length check - HMM will filter garbage
                    int global_dist = prev_gene.distance_from_repeat + c.distance;
                    
                    auto result = ScoreStartNodeWithProfile(
                        c.node, global_dist, c.offset, profile_idx, direction, true);
                    
                    if (result.normalized_score >= (result.is_complete ? params_.MIN_NORMALIZED_SCORE : 0.5)) {
                        result.gene_family = family;
                        result.is_putative = true;
                        if (result.normalized_score > local_best_p2.normalized_score) {
                            local_best_p2 = result;
                        }
                        if (result.normalized_score >= tuning_.phase2_high_confidence) {
                            found_high_confidence_p2.store(true, std::memory_order_relaxed);
                            break;
                        }
                    }
                }

                if (local_best_p2.is_complete || local_best_p2.normalized_score >= PartialScoreThreshold(local_best_p2.profile_index)) {
                    std::lock_guard<std::mutex> lock(p2_mutex);
                    if (local_best_p2.normalized_score > best_p2.normalized_score) {
                        best_p2 = local_best_p2;
                    }
                }
            }
            
            if ((best_p2.is_complete || best_p2.normalized_score >= PartialScoreThreshold(best_p2.profile_index)) && best_p2.normalized_score > next_gene.normalized_score) {
                next_gene = best_p2;
            }
        }

        // ---------------------------------------------------------------------
        // Add gene if found, otherwise stop
        // ---------------------------------------------------------------------
        if (next_gene.normalized_score < (next_gene.is_complete ? params_.MIN_NORMALIZED_SCORE : PartialScoreThreshold(next_gene.profile_index))) {
            std::cout << "  DEBUG: No gene found (is_complete=" << next_gene.is_complete 
                      << " score=" << next_gene.normalized_score << ")" << std::endl;
            cassette.stop_reason_code = "NO_NEXT_GENE";
            break;
        }

        // Post-output length validation: verify amino acid length vs profile bounds
        if (next_gene.profile_index < profiles_.size()) {
            int aa_len = static_cast<int>(next_gene.amino_acids.size());
            const auto& prof = profiles_[next_gene.profile_index];
            if (aa_len < prof.min_aa || aa_len > prof.max_aa) {
                std::cout << "  WARNING: Gene " << next_gene.gene_family
                          << " AA length " << aa_len << " outside profile bounds ["
                          << prof.min_aa << "," << prof.max_aa << "] - skipped" << std::endl;
                cassette.stop_reason_code = "AA_LENGTH_FAIL";
                break;
            }
        }

        // Verify intergenic gap is biologically plausible
        // delta = distance from prev_gene.start_node to next_gene.start_node
        // For UPSTREAM: delta = gap + next_gene.gene_length
        // For DOWNSTREAM: delta = gap (gene extends away from prev)
        int delta = next_gene.distance_from_repeat - prev_gene.distance_from_repeat;
        int actual_gap;
        if (direction == SearchDirection::UPSTREAM) {
            actual_gap = delta - next_gene.gene_length;
        } else {
            // DOWNSTREAM: prev_gene extends toward next_gene
            // delta = prev_gene.gene_length + gap
            actual_gap = delta - prev_gene.gene_length;
        }
        
        if (actual_gap < intergenic_min || actual_gap > intergenic_max) {
            std::cout << "  DEBUG: Gap out of range! delta=" << delta << " actual_gap=" << actual_gap
                      << " allowed=[" << intergenic_min << "," << intergenic_max << "]" << std::endl;
            cassette.stop_reason_code = "GAP_FAIL";
            break;
        }


        cassette.genes.push_back(next_gene);
        cassette.total_distance_bp = next_gene.distance_from_repeat + next_gene.gene_length;
        detected_families.insert(next_gene.gene_family);
        
        // Narrow candidate types based on detected gene family
        // BUT: keep at least one type to avoid losing all context
        if (!candidate_type_indices.empty() && !next_gene.gene_family.empty()) {
            auto narrowed = NarrowCandidateTypesByFamily(
                candidate_type_indices, next_gene.gene_family);
            // Only narrow if we still have candidates; otherwise keep searching broadly
            if (!narrowed.empty()) {
                candidate_type_indices = narrowed;
            }
            // If narrowed is empty, keep the old candidate_type_indices
            // This handles unexpected genes without killing the search
        }
        
        // Update detected_type if narrowed to a single candidate and the rule is satisfied.
        if (candidate_type_indices.size() == 1) {
            const auto& rule = type_rules_[candidate_type_indices[0]];
            auto [mandatory_matches, total_matches] = ScoreTypeMatch(candidate_type_indices[0], detected_families);
            if (mandatory_matches >= rule.min_mandatory && total_matches >= rule.min_genes) {
                cassette.detected_type = rule.type_class;
            }
        } else if (candidate_type_indices.empty() && cassette.detected_type.empty()) {
            cassette.detected_type = "Unknown";
        }
    }

    // Final type determination - score all candidate types and pick best.
    // Score = mandatory_matches * 1000 + total_matches * 10 + log(db_count+1)
    // The log-prior term breaks ties in favour of more prevalent types in CRISPRCasDB.
    if (cassette.detected_type.empty()) {
        auto score_candidates = [&](const std::vector<size_t>& indices) -> bool {
            double best_score = -1.0;
            size_t best_idx   = 0;
            for (size_t type_idx : indices) {
                auto [mandatory_matches, total_matches] = ScoreTypeMatch(type_idx, detected_families);
                const auto& rule = type_rules_[type_idx];
                if (mandatory_matches >= rule.min_mandatory && total_matches >= rule.min_genes) {
                    double sc = mandatory_matches * 1000.0
                              + total_matches     *   10.0
                              + std::log(static_cast<double>(rule.db_count) + 1.0);
                    if (sc > best_score) {
                        best_score = sc;
                        best_idx   = type_idx;
                    }
                }
            }
            if (best_score >= 0.0) {
                cassette.detected_type = type_rules_[best_idx].type_class;
                return true;
            }
            return false;
        };

        if (!candidate_type_indices.empty()) {
            if (!score_candidates(candidate_type_indices)) {
                cassette.detected_type = "Unknown";
            }
        } else {
            auto all_candidates = FindCandidateTypes(detected_families);
            if (!score_candidates(all_candidates)) {
                cassette.detected_type = "Unknown";
            }
        }
    }

    cassette.reached_limit = (cassette.total_distance_bp >= max_cassette_bp);
    if (cassette.reached_limit || (cassette.genes.size() > 1 && cassette.stop_reason_code.empty())) {
        cassette.stop_reason_code = "LIMIT_REACHED";
    }
    
    // CRITICAL: For upstream cassettes, reverse gene order to maintain biological operon direction
    // Upstream search finds closest gene first (near repeat), but biological operon starts farthest
    // Example: gene3 <-- gene2 <-- gene1 <-- REPEAT
    //   Discovery order: gene1, gene2, gene3 (closest to farthest)
    //   Biological order: gene3, gene2, gene1 (operon start to end)
    if (direction == SearchDirection::UPSTREAM && !cassette.genes.empty()) {
        std::reverse(cassette.genes.begin(), cassette.genes.end());
    }
    
    return cassette;
}

// Detect cassettes on BOTH sides of the repeat
std::vector<CasCassette> CasWorkflow::DetectAllCassettes(uint64_t repeat_node) {
    std::cout << __func__ << "\n";
    std::vector<CasCassette> cassettes;

    // Run directions sequentially so each DetectCassette call can use full OpenMP parallelism
    // internally (nested parallel regions are typically disabled and otherwise serialize hot loops).
    CasCassette upstream = DetectCassette(repeat_node, SearchDirection::UPSTREAM);
    CasCassette downstream = DetectCassette(repeat_node, SearchDirection::DOWNSTREAM);

    if (!upstream.genes.empty()) {
        cassettes.push_back(upstream);
    }
    if (!downstream.genes.empty()) {
        cassettes.push_back(downstream);
    }
    
    return cassettes;
}

// Detect cassettes for all repeat nodes from cycles_map
// After post-processing: uses repeat path first/last nodes as BFS origins,
// and marks spacer nodes as blocked so BFS does not traverse into spacer regions.
std::unordered_map<uint64_t, std::vector<CasCassette>> CasWorkflow::DetectAllCassettes(
        const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles_map) {
    std::cout << __func__ << " (cycles_map overload, " << cycles_map.size() << " repeat nodes)\n";

    // Build global bitvector of all cycle nodes (repeat + spacer) so traversals skip them
    cycle_nodes_.assign(sdbg_.size(), false);
    for (const auto& [start_node, cycle_list] : cycles_map) {
        for (const auto& cycle : cycle_list) {
            for (uint64_t node_id : cycle) {
                cycle_nodes_[node_id] = true;
            }
        }
    }

    std::unordered_map<uint64_t, std::vector<CasCassette>> results;

    // Repeat frequency threshold (same as CRISPRPostProcessor)
    constexpr double REPEAT_THRESHOLD = 0.85;

    for (const auto& [repeat_node, cycle_list] : cycles_map) {
        if (cycle_list.size() < 2) {
            // Need at least 2 cycles to compute frequencies
            auto cassettes = DetectAllCassettes(repeat_node);
            if (!cassettes.empty()) results[repeat_node] = std::move(cassettes);
            continue;
        }

        // ─── Step 1: Node frequency f(v) = |{cycles containing v}| / |cycles| ───
        std::unordered_map<uint64_t, int> nodeCount;
        for (const auto& cycle : cycle_list) {
            std::unordered_set<uint64_t> unique(cycle.begin(), cycle.end());
            for (uint64_t v : unique) {
                nodeCount[v]++;
            }
        }

        // ─── Step 2: Repeat set R = { v | f(v) >= tau } ───
        const int threshold = static_cast<int>(
            std::ceil(REPEAT_THRESHOLD * static_cast<double>(cycle_list.size())));
        std::unordered_set<uint64_t> repeatSet;
        for (const auto& [node, count] : nodeCount) {
            if (count >= threshold) {
                repeatSet.insert(node);
            }
        }

        if (repeatSet.empty()) {
            auto cassettes = DetectAllCassettes(repeat_node);
            if (!cassettes.empty()) results[repeat_node] = std::move(cassettes);
            continue;
        }

        // ─── Step 3: Find r0 = node in R with no repeat-internal predecessor ───
        std::unordered_set<uint64_t> successorsOfR;
        for (uint64_t u : repeatSet) {
            uint64_t outgoings[4];
            int numOut = sdbg_.OutgoingEdges(u, outgoings);
            for (int i = 0; i < numOut; ++i) {
                if (repeatSet.count(outgoings[i])) {
                    successorsOfR.insert(outgoings[i]);
                }
            }
        }

        uint64_t r0 = repeat_node;  // fallback
        bool foundStart = false;
        for (uint64_t v : repeatSet) {
            if (successorsOfR.find(v) == successorsOfR.end()) {
                r0 = v;
                foundStart = true;
                break;
            }
        }
        if (!foundStart) {
            r0 = repeatSet.count(repeat_node) ? repeat_node : *repeatSet.begin();
        }

        // ─── Step 4: Build repeat path from longest cycle containing r0 ───
        size_t bestIdx = SIZE_MAX, bestLen = 0;
        for (size_t ci = 0; ci < cycle_list.size(); ++ci) {
            if (std::find(cycle_list[ci].begin(), cycle_list[ci].end(), r0) != cycle_list[ci].end()) {
                if (cycle_list[ci].size() > bestLen) {
                    bestLen = cycle_list[ci].size();
                    bestIdx = ci;
                }
            }
        }

        if (bestIdx == SIZE_MAX) {
            // r0 not found in any cycle — fallback to old single-node behavior
            auto cassettes = DetectAllCassettes(repeat_node);
            if (!cassettes.empty()) results[repeat_node] = std::move(cassettes);
            continue;
        }

        const auto& refCycle = cycle_list[bestIdx];
        size_t r0Pos = static_cast<size_t>(std::distance(refCycle.begin(),
            std::find(refCycle.begin(), refCycle.end(), r0)));

        // Rotate cycle to start at r0
        std::vector<uint64_t> rotatedRef;
        rotatedRef.reserve(refCycle.size());
        for (size_t i = 0; i < refCycle.size(); ++i)
            rotatedRef.push_back(refCycle[(r0Pos + i) % refCycle.size()]);

        // Take leading contiguous repeat-node run as repeat path
        std::vector<uint64_t> repeatPath;
        for (size_t i = 0; i < rotatedRef.size(); ++i) {
            if (repeatSet.count(rotatedRef[i]))
                repeatPath.push_back(rotatedRef[i]);
            else
                break;
        }

        if (repeatPath.empty()) {
            auto cassettes = DetectAllCassettes(repeat_node);
            if (!cassettes.empty()) results[repeat_node] = std::move(cassettes);
            continue;
        }

        // ─── FIRST node for UPSTREAM, LAST node for DOWNSTREAM ───
        uint64_t first_node = repeatPath.front();
        uint64_t last_node  = repeatPath.back();

        std::cout << "  repeat_node=" << repeat_node
                  << " repeat_path_nodes=" << repeatPath.size()
                  << " first_node=" << first_node
                  << " last_node=" << last_node << std::endl;

        // ─── Mark spacer nodes as blocked ───
        // Spacer nodes = non-repeat nodes in every cycle for this group.
        // They are already in cycle_nodes_ from the initial pass above,
        // so BFS will not traverse into spacer regions.

        // ─── Search: UPSTREAM from first repeat node, DOWNSTREAM from last ───
        std::vector<CasCassette> cassettes;

        CasCassette upstream = DetectCassette(first_node, SearchDirection::UPSTREAM);
        upstream.repeat_node = repeat_node;  // preserve original key for result mapping

        CasCassette downstream = DetectCassette(last_node, SearchDirection::DOWNSTREAM);
        downstream.repeat_node = repeat_node;

        if (!upstream.genes.empty()) {
            cassettes.push_back(upstream);
        }
        if (!downstream.genes.empty()) {
            cassettes.push_back(downstream);
        }

        if (!cassettes.empty()) {
            results[repeat_node] = std::move(cassettes);
        }
    }

    return results;
}

// Detect cassettes using POST-PROCESSED filtered arrays.
// Uses the already-computed repeat paths and spacer nodes from CRISPRPostProcessor.
// - repeat_path.front() → UPSTREAM search origin
// - repeat_path.back()  → DOWNSTREAM search origin
// - All repeat + spacer nodes are blocked in cycle_nodes_
std::vector<std::pair<std::string, std::vector<CasCassette>>> CasWorkflow::DetectAllCassettesFromFiltered(
        const std::vector<CRISPRPostProcessor::FilteredArray>& filtered_arrays) {
    std::cout << __func__ << " (" << filtered_arrays.size() << " filtered arrays)\n";

    // Build global bitvector: block ALL repeat + spacer nodes from filtered arrays
    cycle_nodes_.assign(sdbg_.size(), false);
    for (const auto& fa : filtered_arrays) {
        for (uint64_t node_id : fa.repeat_path) {
            cycle_nodes_[node_id] = true;
        }
        for (const auto& sp_path : fa.spacer_node_paths) {
            for (uint64_t node_id : sp_path) {
                cycle_nodes_[node_id] = true;
            }
        }
    }

    std::vector<std::pair<std::string, std::vector<CasCassette>>> results;

    for (const auto& fa : filtered_arrays) {
        if (fa.repeat_path.empty()) continue;

        uint64_t first_node = fa.repeat_path.front();  // UPSTREAM origin
        uint64_t last_node  = fa.repeat_path.back();   // DOWNSTREAM origin

        std::cout << "  repeat=\"" << fa.repeat_sequence.substr(0, 30)
                  << (fa.repeat_sequence.size() > 30 ? "..." : "")
                  << "\" repeat_path_nodes=" << fa.repeat_path.size()
                  << " spacer_paths=" << fa.spacer_node_paths.size()
                  << " first_node=" << first_node
                  << " last_node=" << last_node << std::endl;

        std::vector<CasCassette> cassettes;

        CasCassette upstream = DetectCassette(first_node, SearchDirection::UPSTREAM);
        upstream.repeat_node = first_node;

        CasCassette downstream = DetectCassette(last_node, SearchDirection::DOWNSTREAM);
        downstream.repeat_node = last_node;

        if (!upstream.genes.empty()) {
            cassettes.push_back(upstream);
        }
        if (!downstream.genes.empty()) {
            cassettes.push_back(downstream);
        }

        if (!cassettes.empty()) {
            results.push_back({fa.repeat_sequence, std::move(cassettes)});
        }
    }

    return results;
}
