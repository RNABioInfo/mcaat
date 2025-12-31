# Changelog

## [Unreleased] - 0.5.1
- Fixed crash in spacer ordering when no reads are found (guard added).
- Improved parallel scaling:
  - Replaced frequent `#pragma omp critical` usage with per-thread buffers and serial merges.
  - Introduced a lock-free visited bitmap (atomic 64-bit words) to remove synchronization hot-spots.
- Reduced allocator contention and reused per-thread containers to lower memory churn under heavy parallelism.
- Build fixes: added missing includes and small portability fixes so CMake build succeeds on target platforms.
- Build: CMake configure and full build completed successfully (targets `mcaat` and `runTests` built).


*Notes:* these changes focus on improving scalability and preventing serialization on large, memory-bound graphs.
