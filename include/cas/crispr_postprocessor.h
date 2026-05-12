#ifndef CRISPR_POSTPROCESSOR_H
#define CRISPR_POSTPROCESSOR_H

#include <string>
#include <vector>
#include <cstdint>

/**
 * @brief Bridging types between the CRISPR array detection pipeline
 *        (PostProcessor) and the Cas gene detection (CasWorkflow).
 *
 * CRISPRPostProcessor::FilteredArray holds all data produced by the
 * post-processing step that CasWorkflow needs to run cassette detection:
 *   - repeat_path          : ordered SDBG node IDs spanning the repeat
 *   - spacer_node_paths    : per-spacer ordered SDBG node IDs
 *   - repeat_sequence      : consensus repeat nucleotide sequence
 */
namespace CRISPRPostProcessor {

struct FilteredArray {
    std::string repeat_sequence;                          // Consensus repeat
    std::vector<uint64_t> repeat_path;                   // SDBG nodes covering the repeat
    std::vector<std::vector<uint64_t>> spacer_node_paths; // SDBG nodes per spacer
};

}  // namespace CRISPRPostProcessor

#endif  // CRISPR_POSTPROCESSOR_H
