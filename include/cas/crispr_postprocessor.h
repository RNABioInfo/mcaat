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
 *   - repeat_sequence      : supported observed repeat nucleotide sequence
 */
namespace CRISPRPostProcessor {

struct FilteredArray {
    std::string repeat_sequence;                          // Supported observed repeat (legacy files may contain a consensus)
    std::vector<uint64_t> repeat_path;                   // Valid SDBG edge handles whose source k-mers cover the repeat
    std::vector<std::vector<uint64_t>> spacer_node_paths; // Valid SDBG edge handles for pure spacer k-mers
};

}  // namespace CRISPRPostProcessor

#endif  // CRISPR_POSTPROCESSOR_H
