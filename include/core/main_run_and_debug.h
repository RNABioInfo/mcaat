#ifndef INCLUDE_MAIN_RUN_AND_DEBUG_H_
#define INCLUDE_MAIN_RUN_AND_DEBUG_H_

#include <tuple>
#include <vector>
#include <chrono>

#include "core/settings.h"
#include "sdbg/sdbg.h"
#include "io/reads.h"
#include "crispr_array/evaluation.h"
#include "crispr_array/spacer_ordering.h"
#include "core/tmp_utils.h"

vector<vector<uint64_t>>
run_and_debug_finding_of_relevant_reads(
    const vector<vector<uint64_t>>& cycles,
    const Settings& settings,
    const SDBG& sdbg
);

vector<tuple<string, string, vector<string>, float, float>>
run_and_debug_spacer_ordering(
    const vector<vector<uint64_t>>& reads,
    SDBG& sdbg,
    const vector<vector<uint64_t>>& cycles
);

void run_and_debug_benchmark_results(
    const Settings& settings,
    const vector<tuple<string, string, vector<string>, float, float>>& found_systems
);

void run_and_debug_results(
    const vector<tuple<string, string, vector<string>, float, float>>& found_systems
);

#endif
