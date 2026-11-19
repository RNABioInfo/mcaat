# Changelog

## Protospacer workflow v10 / seeded v7
- Source-array crossing is off by default; it followed the array into a neighbouring spacer's phage.
- New CRISPR repeat barrier (`--crispr-repeat-barrier true`): flanks stop at any reported repeat.
- New `scripts/evaluate_target_region.py`: counts a match only when the reference region contains the target.
- New `tests/check_array_chimera.py` and `scripts/run_evidence19_stages.sh`. See `docs/protospacer-target-region.md`.

## [Unreleased] - 0.5.1
- Fixed crash in spacer ordering when no reads are found (guard added).
- Improved parallel scaling:
  - Replaced frequent `#pragma omp critical` usage with per-thread buffers and serial merges.
  - Introduced a lock-free visited bitmap (atomic 64-bit words) to remove synchronization hot-spots.
- Reduced allocator contention and reused per-thread containers to lower memory churn under heavy parallelism.
- Build fixes: added missing includes and small portability fixes so CMake build succeeds on target platforms.
- Build: CMake configure and full build completed successfully (targets `mcaat` and `runTests` built).


*Notes:* these changes focus on improving scalability and preventing serialization on large, memory-bound graphs.
