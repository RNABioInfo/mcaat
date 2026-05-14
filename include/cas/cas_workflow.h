#ifndef CAS_WORKFLOW_H
#define CAS_WORKFLOW_H

#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cstdint>
#include <cmath>
#include "sdbg/sdbg.h"
#include "profile.h"
#include "buckets.h"

// Direction for cassette detection relative to CRISPR repeat
enum class SearchDirection {
    DOWNSTREAM,  // CAS genes after repeat (original behavior)
    UPSTREAM     // CAS genes before repeat (reverse search)
};

struct CasWorkflowParams {
    int FIRST_GENE_MIN_DIST = 50;
    int FIRST_GENE_MAX_DIST = 2500;
    int MAX_START_CANDIDATES = 50000;
    int BEAM_WIDTH = 100;
    double MIN_NORMALIZED_SCORE = 0.2;  // bits per HMM position; rejects random matches
    bool allow_exploratory = true;  // If false, stop when rule-guided search fails (no Phase 2)
};

// Extract base gene family from profile name for deduplication
// New naming convention: "cas3_I_1" -> "cas3", "cas8a1_I-A_1" -> "cas8a1"
// Old naming convention: "Cas3HD_0_I_1" -> "cas3", "Cas1_0_I-A_1" -> "cas1" (normalized)
inline std::string ExtractGeneFamily(const std::string& profile_name) {
    // Remove .hmm extension if present
    std::string name = profile_name;
    size_t ext_pos = name.rfind(".hmm");
    if (ext_pos != std::string::npos) {
        name = name.substr(0, ext_pos);
    }
    
    // Normalize to lowercase for consistency
    std::string lower_name;
    lower_name.reserve(name.size());
    for (char c : name) {
        lower_name += std::tolower(c);
    }
    
    // Extract first part before '_'
    size_t underscore = lower_name.find('_');
    std::string base = (underscore != std::string::npos) ? lower_name.substr(0, underscore) : lower_name;
    
    // Handle domain profiles before the generic CAS parser.
    if (base.rfind("cas3hd", 0) == 0) {
        return "cas3";
    }
    if (base.rfind("cas10hd", 0) == 0) {
        return "cas10";
    }
    if (base.rfind("cas10-like", 0) == 0) {
        return "cas10-like";
    }

    // Keep cas prefix + digits + optional letter suffix (a/b/c/d), strip other domain suffixes
    if (base.size() >= 3 && base.substr(0, 3) == "cas") {
        size_t i = 3;
        // Include digits
        while (i < base.size() && std::isdigit(base[i])) {
            ++i;
        }
        // Check for letter suffixes like 'a', 'b', 'c', 'd' (e.g., cas8a, cas8b)
        if (i < base.size() && std::isalpha(base[i]) && std::islower(base[i])) {
            ++i;  // include the lowercase letter
            // Include trailing digit if present (e.g., cas8a1, cas8b12)
            while (i < base.size() && std::isdigit(base[i])) {
                ++i;
            }
        }
        return base.substr(0, i);
    }
    
    // For non-Cas genes (csa5, csax, csm3, etc.), extract up to first underscore
    return base;
}

struct DetectedCasGene {
    std::string profile_name;
    std::string gene_family;     // Base family for dedup: Cas1, Cas3, Cas8a1...
    size_t profile_index = SIZE_MAX;  // Index into profiles_ array
    uint64_t start_node = 0;
    uint64_t end_node = 0;
    int distance_from_repeat = 0;
    int gene_length = 0;
    double normalized_score = -1e9;
    double bit_score = 0.0;
    double evalue = -1.0;  // E-value from Gumbel EVD (-1 = not computed)
    std::string amino_acids;
    std::vector<uint64_t> node_path;
    bool is_complete = false;
    bool is_putative = false;    // True if found in exploratory phase
};

struct CasCassette {
    uint64_t repeat_node = 0;
    std::vector<DetectedCasGene> genes;
    int total_distance_bp = 0;
    bool reached_limit = false;
    std::string stop_reason_code = "";  // NO_NEXT_GENE | GAP_FAIL | LIMIT_REACHED
    std::string detected_type;   // e.g., "I-A", "I-C", "II-A"
    SearchDirection direction = SearchDirection::DOWNSTREAM;  // Which side of repeat
};

struct StartCodonCandidate {
    uint64_t node;
    int distance;
    int offset;
    int orf_length = 0;  // Pre-computed ORF length for filtering
};

// Result of testing one CRISPR type hypothesis
struct TypeHypothesisResult {
    std::string type_class;
    std::vector<DetectedCasGene> genes;
    double total_score = 0.0;
    int mandatory_found = 0;
    int total_genes = 0;
    int db_count = 1;  // Propagated from CasTypeRule for log-prior scoring
    
    double CassetteScore() const {
        // Matches live scoring in score_candidates lambda:
        // mandatory_found dominates; total_genes breaks ties; log-prior favours common types
        return mandatory_found * 1000.0
             + total_genes    *   10.0
             + std::log(static_cast<double>(db_count) + 1.0);
    }
};

// Rule loaded from _rules.csv (based on XML definitions from MacSyFinder CasFinder)
// CSV format: type_class,mandatory_profiles,accessory_profiles,min_mandatory,min_genes
// mandatory = type-defining genes (at least min_mandatory must be found)
// accessory = optional genes including adaptation (cas1/2/4) and others
struct CasTypeRule {
    std::string type_class;                       // e.g., "I-A", "II-C", "V-A"
    std::vector<std::string> mandatory_genes;     // Type-defining genes (was interference_genes)
    std::vector<std::string> accessory_genes;     // Optional genes (includes adaptation + others)
    
    // Gene families extracted from profile names for family-based matching
    // Each entry is a set of alternative families (e.g., {"cas8a"} or {"cas5", "csc1"})
    std::vector<std::set<std::string>> mandatory_families;
    std::vector<std::set<std::string>> accessory_families;
    
    int min_mandatory = 1;                        // Minimum mandatory genes required
    int min_genes = 1;                            // Minimum total genes required
    int max_genes = 20;                           // Maximum total genes allowed
    int db_count  = 1;                            // Verified instances in CRISPRCasDB (log-prior)
    
    // Check if a detected gene family matches any of the mandatory genes
    bool MatchesMandatoryFamily(const std::string& family) const {
        for (const auto& alt_set : mandatory_families) {
            if (alt_set.count(family) > 0) return true;
        }
        return false;
    }
    
    // Check if a detected gene family matches any of the accessory genes
    bool MatchesAccessoryFamily(const std::string& family) const {
        for (const auto& alt_set : accessory_families) {
            if (alt_set.count(family) > 0) return true;
        }
        return false;
    }
    
    // Backwards compatibility aliases
    bool MatchesInterferenceFamily(const std::string& family) const {
        return MatchesMandatoryFamily(family);
    }
    bool MatchesAdaptationFamily(const std::string& family) const {
        return MatchesAccessoryFamily(family);
    }
    
    // Count how many mandatory gene slots are satisfied by a set of detected families
    int CountMatchedMandatory(const std::set<std::string>& detected_families) const {
        int count = 0;
        for (const auto& alt_set : mandatory_families) {
            for (const auto& fam : alt_set) {
                if (detected_families.count(fam) > 0) {
                    ++count;
                    break;  // This slot is satisfied, move to next
                }
            }
        }
        return count;
    }
    
    // Count how many total gene slots are satisfied
    int CountMatchedTotal(const std::set<std::string>& detected_families) const {
        int count = CountMatchedMandatory(detected_families);
        for (const auto& alt_set : accessory_families) {
            for (const auto& fam : alt_set) {
                if (detected_families.count(fam) > 0) {
                    ++count;
                    break;
                }
            }
        }
        return count;
    }
    
    // Backwards compatibility alias
    int CountMatchedInterference(const std::set<std::string>& detected_families) const {
        return CountMatchedMandatory(detected_families);
    }
    
    // Find position of a gene family in the rule's gene order (mandatory + accessory)
    // Returns -1 if not found
    // Gene order: [mandatory_genes..., accessory_genes...]
    int FindFamilyPosition(const std::string& family) const {
        // Check mandatory first
        for (size_t i = 0; i < mandatory_families.size(); ++i) {
            if (mandatory_families[i].count(family) > 0) {
                return static_cast<int>(i);
            }
        }
        // Then accessory
        for (size_t i = 0; i < accessory_families.size(); ++i) {
            if (accessory_families[i].count(family) > 0) {
                return static_cast<int>(mandatory_families.size() + i);
            }
        }
        return -1;
    }
    
    // Get families that come BEFORE position (for UPSTREAM chaining)
    // Returns set of families at positions 0 to pos-1
    std::set<std::string> GetFamiliesBeforePosition(int pos) const {
        std::set<std::string> result;
        int total_mandatory = static_cast<int>(mandatory_families.size());
        
        for (int i = 0; i < pos; ++i) {
            if (i < total_mandatory) {
                for (const auto& fam : mandatory_families[i]) {
                    result.insert(fam);
                }
            } else {
                int acc_idx = i - total_mandatory;
                if (acc_idx < static_cast<int>(accessory_families.size())) {
                    for (const auto& fam : accessory_families[acc_idx]) {
                        result.insert(fam);
                    }
                }
            }
        }
        return result;
    }
    
    // Get families that come AFTER position (for DOWNSTREAM chaining)
    // Returns set of families at positions pos+1 to end
    std::set<std::string> GetFamiliesAfterPosition(int pos) const {
        std::set<std::string> result;
        int total_mandatory = static_cast<int>(mandatory_families.size());
        int total = total_mandatory + static_cast<int>(accessory_families.size());
        
        for (int i = pos + 1; i < total; ++i) {
            if (i < total_mandatory) {
                for (const auto& fam : mandatory_families[i]) {
                    result.insert(fam);
                }
            } else {
                int acc_idx = i - total_mandatory;
                if (acc_idx < static_cast<int>(accessory_families.size())) {
                    for (const auto& fam : accessory_families[acc_idx]) {
                        result.insert(fam);
                    }
                }
            }
        }
        return result;
    }
};

#include "crispr_postprocessor.h"

class CasWorkflow {
public:
    CasWorkflow(SDBG& sdbg, const std::string& profiles_dir, const std::string& rules_csv_path);
    void ApplySensitivityMode(bool enabled = true);
    
    // Detect cassettes on BOTH sides of repeat (upstream + downstream)
    std::vector<CasCassette> DetectAllCassettes(uint64_t repeat_node);
    
    // Detect cassettes for all repeat nodes from cycles_map (raw cycles — legacy)
    std::unordered_map<uint64_t, std::vector<CasCassette>> DetectAllCassettes(
        const std::unordered_map<uint64_t, std::vector<std::vector<uint64_t>>>& cycles_map);
    
    // Detect cassettes using POST-PROCESSED filtered arrays
    // Uses repeat_path.front() for UPSTREAM, repeat_path.back() for DOWNSTREAM
    // All repeat + spacer nodes are blocked so BFS cannot enter those regions
    std::vector<std::pair<std::string, std::vector<CasCassette>>> DetectAllCassettesFromFiltered(
        const std::vector<CRISPRPostProcessor::FilteredArray>& filtered_arrays);
    
    // Detect cassette in a specific direction (internal method, also exposed for testing)
    CasCassette DetectCassette(uint64_t repeat_node, SearchDirection direction);
    
    // Legacy method - detects only downstream (for backwards compatibility)
    std::vector<DetectedCasGene> DetectCasGenes(uint64_t repeat_node);

    void SetParams(const CasWorkflowParams& params) { params_ = params; }
    const CasWorkflowParams& GetParams() const { return params_; }
    const std::vector<HMMProfiles::ProfileSize>& GetProfiles() const { return profiles_; }

private:
    struct WorkflowTuning {
        int first_gene_search_max_bp = 8000;
        int bfs_max_candidates = 10000;
        int max_cassette_bp = 30000;
        int intergenic_min_bp = -24;
        int intergenic_max_bp = 2000;
        int max_genes = 10;
        int first_gene_bins = 20;
        int profile_window_count = 50;
        double shallow_threshold = 0.05;
        double phase1_high_confidence = 1.2;
        double phase2_high_confidence = 0.85;
    };

    SDBG& sdbg_;
    std::string profiles_dir_;
    CasWorkflowParams params_;
    WorkflowTuning tuning_;
    bool sensitivity_mode_ = false;
    std::vector<HMMProfiles::ProfileSize> profiles_;
    
    // Rules loaded from _rules.csv
    std::vector<CasTypeRule> type_rules_;
    
    // First-gene profile INDICES (not filenames) - all profiles that can start a cassette
    std::vector<size_t> first_gene_profile_indices_;
    
    // First-gene REPRESENTATIVES: one (longest) profile per family for fast initial scan
    // Format: (profile_index, family_name)
    std::vector<std::pair<size_t, std::string>> first_gene_representatives_;
    
    // Map: profile index -> which type rules use this as first gene
    std::map<size_t, std::vector<std::string>> first_gene_to_types_;
    
    // Profiles storage - loaded lazily with mutex protection
    mutable std::vector<Profile> loaded_profiles_;
    mutable std::vector<char> profile_loaded_;  // Track which are loaded (char, not bool — vector<bool> is bit-packed and not thread-safe)
    
    const Profile* GetProfile(size_t profile_index) const;
    
    // Maps gene family -> list of profile indices for that family
    std::map<std::string, std::vector<size_t>> family_to_profiles_;
    
    // HMM score cache: (start_node, profile_index, offset) -> DetectedCasGene
    // Eliminates redundant HMM computations across overlapping search windows
    mutable std::map<std::tuple<uint64_t, size_t, int>, DetectedCasGene> score_cache_;
    
    // Bitvector of cycle/spacer nodes to skip during BFS traversals
    std::vector<bool> cycle_nodes_;
    
    std::string GetNodeSequence(uint64_t node_id);
    bool HasStartCodon(uint64_t node_id);
    int GetStartCodonOffset(uint64_t node_id);
    bool HasStopCodon(uint64_t node_id);
    int GetStopCodonOffset(uint64_t node_id);
    int EstimateORFLength(uint64_t start_node, int start_offset, int max_search);
    int EstimateORFLengthReverse(uint64_t stop_node, int stop_offset, int max_search);
    bool IsProfileCompatibleWithDistance(const HMMProfiles::ProfileSize& profile, int distance_from_repeat, SearchDirection direction);

    // Returns the acceptable incomplete-alignment score threshold for a profile.
    // Large profiles (LENG > 500) get a lower threshold: banded Viterbi often
    // drifts on multi-domain proteins and fails to consume the full HMM.
    double PartialScoreThreshold(size_t profile_index) const {
        if (profile_index < profiles_.size() && profiles_[profile_index].leng > 500)
            return 0.3;
        return 0.5;
    }
    
    std::vector<StartCodonCandidate> FindStartCodonCandidates(
        uint64_t repeat_node, int min_dist, int max_dist, int max_candidates,
        SearchDirection direction);
    
    std::vector<StartCodonCandidate> FindStartCodonCandidatesFromNode(
        uint64_t origin_node, int min_dist, int max_dist, int max_candidates,
        SearchDirection direction);
    
    DetectedCasGene ScoreStartNodeWithProfile(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        size_t profile_index, SearchDirection direction,
        bool skip_distance_check = false);
    
    // PROGRESSIVE DEPTH: Quick shallow scan to reject non-promising candidates
    double ScoreStartNodeShallow(
        uint64_t start_node, int start_codon_offset,
        size_t profile_index, SearchDirection direction);
    
    DetectedCasGene ScoreStartNodeWithAllProfiles(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        SearchDirection direction);
    
    // Score with only profiles matching expected families (rule-guided phase)
    DetectedCasGene ScoreStartNodeWithExpectedFamilies(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        const std::set<std::string>& expected_families, SearchDirection direction);
    
    // Score with all profiles but exclude already-detected families (exploratory phase)
    DetectedCasGene ScoreStartNodeExploratory(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        const std::set<std::string>& detected_families, SearchDirection direction);
    
    // MULTI-PROFILE SCORING: Score ONE candidate against ALL specified profiles
    // in a single graph traversal. Returns vector of (profile_idx, DetectedCasGene) pairs.
    std::vector<std::pair<size_t, DetectedCasGene>> ScoreStartNodeMultiProfile(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        const std::vector<std::pair<size_t, std::string>>& profile_family_pairs,  // (profile_idx, family)
        SearchDirection direction);
    
    // Build family-to-profiles index on startup
    void BuildFamilyIndex();
    
    // Load rules from _rules.csv
    void LoadRulesFromCSV(const std::string& rules_csv_path);
    
    // Build first-gene profile index from rules
    void BuildFirstGeneIndex();
    
    // Score with ONLY the first-gene profiles (from rules) - UPSTREAM
    DetectedCasGene ScoreStartNodeWithFirstGeneProfiles(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        SearchDirection direction);
    
    // Score with first-gene profiles using shallow pre-filtering - DOWNSTREAM only
    DetectedCasGene ScoreStartNodeWithFirstGeneProfilesDownstream(
        uint64_t start_node, int distance_from_repeat, int start_codon_offset,
        SearchDirection direction);
    
    // Get profile indices for genes in a specific type
    std::vector<size_t> GetProfileIndicesForType(const std::string& type_class) const;
    
    // Score how well a set of detected gene families matches a candidate type
    // Returns: pair<interference_matches, total_matches>
    std::pair<int, int> ScoreTypeMatch(
        size_t type_idx, const std::set<std::string>& detected_families) const;
    
    // Find all candidate types that could match a set of detected gene families
    std::vector<size_t> FindCandidateTypes(
        const std::set<std::string>& detected_families) const;
    
    // Narrow candidate types based on a newly detected gene family
    std::vector<size_t> NarrowCandidateTypesByFamily(
        const std::vector<size_t>& candidate_type_indices,
        const std::string& detected_family) const;
    
    // Type-driven detection: test a single type hypothesis
    TypeHypothesisResult DetectCassetteForType(
        const std::vector<StartCodonCandidate>& candidates,
        const CasTypeRule& rule,
        SearchDirection direction);
    
    // Get all profile indices for a type (mandatory + accessory)
    std::vector<std::pair<size_t, bool>> GetTypeProfiles(const CasTypeRule& rule) const;
};

#endif // CAS_WORKFLOW_H
