#ifndef MCAAT_PROTO_ARRAY_CONTEXT_H
#define MCAAT_PROTO_ARRAY_CONTEXT_H

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace mcaat::proto {

// Both strings are oriented outwards from a seed boundary. The unmatched
// spacer portion supplies the distance to the repeat, not an identity test:
// a related spacer in an array must not become a target merely because its
// unaligned bases differ. Repeat copies may contain substitutions/indels.
struct FlankContext {
    std::string spacer_remainder;
    std::string repeat;
};
enum class ContextState { pending, external, array };

inline ContextState context_state(const std::string& outward,
                                  const std::vector<FlankContext>& models,
                                  double identity) {
    bool possible = false;
    for (const auto& model : models) {
        if (outward.size() < model.spacer_remainder.size()) {
            possible = true;
            continue;
        }
        const auto& repeat = model.repeat;
        size_t edits = static_cast<size_t>(std::floor((1 - identity) * repeat.size() + 1e-9));
        size_t available = outward.size() - model.spacer_remainder.size();
        size_t observed = std::min(available, repeat.size() + edits);
        // Global repeat alignment against successive prefixes of the flank.
        // No free clipping of the repeat, spacer, or start of the flank.
        std::vector<size_t> row(repeat.size() + 1), next(row.size());
        for (size_t j = 0; j < row.size(); ++j) row[j] = j;
        for (size_t i = 1; i <= observed; ++i) {
            next[0] = i;
            char base = outward[model.spacer_remainder.size() + i - 1];
            for (size_t j = 1; j < row.size(); ++j)
                next[j] = std::min({row[j] + 1, next[j - 1] + 1,
                                   row[j - 1] + (base != repeat[j - 1])});
            row.swap(next);
            if (i + edits >= repeat.size() && row.back() <= edits)
                return ContextState::array;
        }
        if (available < repeat.size() + edits &&
            *std::min_element(row.begin(), row.end()) <= edits) possible = true;
    }
    return possible ? ContextState::pending : ContextState::external;
}

inline std::vector<FlankContext> flank_contexts(const std::string& query,
                                               size_t begin, size_t end,
                                               const std::vector<std::string>& repeats,
                                               bool backward) {
    std::vector<FlankContext> result;
    for (auto repeat : repeats) {
        auto remainder = backward ? query.substr(0, begin) : query.substr(end);
        if (backward) {
            std::reverse(remainder.begin(), remainder.end());
            std::reverse(repeat.begin(), repeat.end());
        }
        result.push_back({std::move(remainder), std::move(repeat)});
    }
    return result;
}

} // namespace mcaat::proto
#endif
