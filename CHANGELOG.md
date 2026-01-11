# Changelog

## [Unreleased] - 0.5.2
- **Renamed**: AminoAcidator → CasGeneDetector throughout codebase
- **Fixed**: max_depth = HMM_length × 3 (was hardcoded 1200)
- **Fixed**: start_node required via command line (removed 100k-200k search)
- **Validation**: 704 Cas gene profiles, 2.16% average difference vs HMMER
- **Documentation**: Streamlined all .md files (concise, preserved math/algorithms)

## [0.5.1]
- Fixed crash when no reads found (guard added)
- Parallel scaling improvements:
  - Replaced `#pragma omp critical` with per-thread buffers
  - Lock-free visited bitmap (atomic 64-bit)
- Reduced allocator contention
- Build fixes: missing includes, CMake portability
- Build: Successfully compiled `mcaat` and `runTests`
