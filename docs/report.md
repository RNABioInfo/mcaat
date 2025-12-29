# MCAAT: Cycle Finder — Algorithmic & Optimization Report ✅

**Scope:** This document describes *only* the algorithmic changes and optimizations introduced in the `optimizations` branch for the cycle finder logic. It is organized so you can read step-by-step what changed, why it was done, and the expected impact.

---

## Summary

1. Replaced global critical sections and shared writes with *per-thread buffers* and a single serial merge step to remove contention.  
2. Replaced lock-based or synchronized "visited" bookkeeping with a *lock-free atomic bitset* (1 bit per node).  
3. Reduced allocations and allocator contention by *reusing per-thread pools* (megahit-style) and preallocating where useful.  
4. Applied traversal micro-optimizations (fixed-size arrays, prefetch, branch hints) to reduce per-edge overhead.

---

## Step-by-step changes (algorithmic & optimization focus)

1) Remove global critical sections → Per-thread buffers + serial merge 🔧
   - What changed:
     - Replaced OpenMP `#pragma omp critical` style updates to shared containers with `vector<...>` of per-thread collectors (e.g., `local_chunks` and `local_results`).
     - After the parallel loop completes, a single-threaded loop merges per-thread buffers into the shared map or results container.
   - Files/locations:
     - `CycleFinder::ChunkStartNodes` (collect start nodes into `local_chunks[tid]` then merge).
     - `CycleFinder::FindApproximateCRISPRArrays` (collect per-thread `local_results`, then merge into `this->results`).
   - Why / Benefit:
     - Eliminates high-contention points on hot shared data structures, enabling scaling to higher core counts.
     - Serial merge cost is amortized and avoids expensive locking in hot loops.

2) Lock-free visited bitmap (1 bit per node) 🔒→⚡
   - What changed:
     - Introduced a global `std::vector<uint64_t> s_visited_words` as a bitset (one bit per node).
     - Provided helper inline functions: `InitializeVisitedGlobal(n)`, `IsVisitedGlobal(node)`, and `MarkVisitedGlobal(node)` implemented using GCC/Clang atomic builtins (`__atomic_load_n`, `__atomic_fetch_or`) with `__ATOMIC_RELAXED` ordering.
   - Files/locations:
     - `src/cycle_finder.cpp` (static `s_visited_words` and helpers) and uses in `FindCycle`, `FindCycleUtil`, and background checks.
   - Why / Benefit:
     - Avoids `vector<std::atomic>` pitfalls (copy/resize/copyability) and the overhead of locks around visited updates.
     - One atomic word operation per change (bit flip) is cheap and scales well.
     - Memory is compact (1 bit per node) and predictable for large graphs.
   - Correctness note:
     - Using relaxed atomics is acceptable because bits only transition from 0→1 (monotonic); races among writers do not break correctness, and reads can tolerate transient states.

3) Reduce allocations and reuse per-thread pools (megahit-style) ♻️
   - What changed:
     - Introduced `static thread_local` pools for DLS (`dls_stack_pool` and `dls_visited_pool`) used by `DepthLevelSearch`.
     - Pools are `clear()`ed between uses but retain capacity; small initial reserve is set to avoid repeated small allocations.
   - Files/locations:
     - `CycleFinder::DepthLevelSearch`.
   - Why / Benefit:
     - Avoids heavy allocator contention when many threads create/destroy temporaries frequently.
     - Reduced per-edge latency and improved throughput during parallel graph traversal.

4) Traversal micro-optimizations (branch hints, fixed arrays, prefetch) 🧠
   - What changed:
     - Use of fixed-size neighbor arrays (`uint64_t neighbors[MAX_EDGE_COUNT]`) rather than heap allocations per node.
     - Prefetching neighbor buffers and using `__builtin_expect` branch hints to optimize hot paths.
     - Small loop unrolling where out-degree is small (de Bruijn graph pattern) to reduce loop overhead.
   - Files/locations:
     - `DepthLevelSearch`, `_GetOutgoings`, and `_GetIncomings` helpers.
   - Why / Benefit:
     - Better cache locality and fewer branch mispredictions; straightforward per-edge speedups with little code complexity.

5) Results merging and memory hygiene 🧽
   - What changed:
     - Per-thread `local_results` (maps) are merged serially into `this->results` after each bucket processed.
     - Call `malloc_trim(0)` occasionally after buckets to release heap fragments back to the OS (for long runs with variable memory usage).
   - Files/locations:
     - `FindApproximateCRISPRArrays`.
   - Why / Benefit:
     - Reduces concurrent unordered_map modification (expensive) and helps long-running runs avoid growing memory footprints unnecessarily.

---

## Expected performance and behavior improvements

- Improved scalability with thread counts beyond the earlier observed plateau (~24 cores) because:
  - Contention points are removed or drastically reduced.
  - Allocator pressure is lowered by reusing containers.
  - Atomic operations on compact bitmaps replace heavier locks.
- Memory cost: the visited bitset adds ~1 bit per node (compact) and per-thread buffers increase transient memory usage proportional to thread count but only for selected nodes.

---

## Limitations & future work

- NUMA-aware allocation and memory binding were not implemented yet — this is the natural next step for large multi-socket machines where memory bandwidth dominates.
- Further profiling (perf/VTune) is needed to quantify the exact causes of any remaining scalability bottlenecks (cache-line bouncing, allocator hotspots, or procedural serial sections).

---

## How to validate quickly (recommended)

1. Check out the `optimizations` branch.
2. Build (`cmake .. && make -j`) and run the same workload used before.
3. Compare (a) execution time vs thread count (1, 8, 24, 48, 128), (b) throughput (nodes/sec), and (c) cycles found to ensure no correctness regression.
4. Use `perf top` / `perf record` or `numastat` to verify reduced lock/atomic time and identify remaining hotspots.

---

## Files touched (algorithmic/optimization only)

- `src/cycle_finder.cpp` — main implementation of lock-free visited bitmap, per-thread collectors, DLS pools, traversal micro-optimizations, merging logic.
- `include/cycle_finder.h` — updated helpers and declarations related to visited bookkeeping (if applicable).

---

## TL;DR

- Replaced shared locks with per-thread buffers + serial merges, added a compact lock-free visited bitmap, and reduced allocation churn via per-thread pools. These changes reduce contention and allocator pressure and improve multithreaded scaling while keeping memory usage reasonable for very large graphs.

---

### Expanded Documentation: files, API, usage, testing, and benchmarks 🧾

Below are additional sections documenting the implementation and usage in more depth.

---

## Detailed file-by-file changes 🔍

- `src/cycle_finder.cpp`
  - Implemented the lock-free visited bitmap (`s_visited_words`) and helpers:
    - `InitializeVisitedGlobal(size_t n_nodes)`
    - `bool IsVisitedGlobal(uint32_t node)`
    - `void MarkVisitedGlobal(uint32_t node)`
  - Added per-thread collectors and serial merge logic used in `ChunkStartNodes` and `FindApproximateCRISPRArrays`.
  - Implemented `static thread_local` DLS pools used by `DepthLevelSearch`.
  - Added traversal micro-optimizations: prefetch, fixed-size neighbor buffers, and branch hints.
  - Included optional calls to `malloc_trim(0)` in long-running buckets.

- `include/cycle_finder.h`
  - Exposed new helper declarations (initialization and debugging helpers) and configuration knobs for pool sizes and prefetch tuning.

- `tests/` (recommended additions)
  - `tests/visited_bitmap_stress.cpp`: small test that races many threads marking visited bits to validate monotonic behavior.
  - `tests/benchmark_multi_threaded.sh`: a script to run the cycle finder across multiple thread counts and collect wall-clock results.

---

## API reference (important functions) 📚

- `void InitializeVisitedGlobal(size_t n_nodes)`
  - Pre-allocates and zero-initializes the visited bitset for `n_nodes` nodes.

- `bool IsVisitedGlobal(uint32_t node)`
  - Returns whether a node is marked visited. Uses relaxed atomics; can return false positives transiently but never returns a persistent false negative once a node is marked.

- `void MarkVisitedGlobal(uint32_t node)`
  - Atomically sets the visited bit for `node` using `__atomic_fetch_or` with `__ATOMIC_RELAXED`.

- `std::vector<Cycle> CycleFinder::FindApproximateCRISPRArrays(...)`
  - Main entry point for the current algorithm; now uses per-thread collectors and merges results serially.

- `DepthLevelSearch` (internal)
  - Uses per-thread pools and fixed-size buffers; configured via compile-time constants or optional runtime flags in `cycle_finder.h`.

> Note: The above helpers are designed for monotonic bitset updates (0→1 only); they are not a general-purpose concurrent set implementation.

---

## Build & run (quick guide) ⚙️

1. Checkout the branch:
   - `git checkout optimizations`
2. Build (out-of-source recommended):
   - `mkdir -p build && cd build`
   - `cmake .. && make -j$(nproc)`
3. Run the cycle finder with typical arguments (example):
   - `./bin/cycle_finder --input graphs/huge_graph.bin --threads 24 --out results.json`
4. Recommended environment variables:
   - `MALLOC_ARENA_MAX=4` to limit glibc arenas and reduce allocator noise for some workloads.

---

## Benchmarks & validation plan 📈

- Quick smoke tests:
  - Run the binary with `--verify-only` or with small graphs to confirm identical cycle sets are returned versus the baseline branch.

- Automated benchmark script (recommended): `bench/run_benchmarks.sh`
  - Runs the same workload across a list of thread counts (e.g., 1, 4, 8, 24, 48) and records wall clock and cycles found in a CSV.
  - Produces a simple gnuplot/matplotlib script to visualize scaling.

- Profiling:
  - Use `perf record -g` and `perf report` to verify decreased time in atomic/lock hotspots and reduced allocator overhead.
  - On NUMA machines, use `numactl --show` and `numastat` to diagnose memory binding issues.

---

## Testing suggestions ✅

- Unit tests:
  - Add a deterministic `visited_bitmap_unit_test` that runs many threads marking bits and then verifies the expected counts.

- Integration tests:
  - A small-run integration test that compares cycles discovered by the optimized implementation against the canonical reference implementation for multiple small graphs.

- Stress tests:
  - Large graph runs with tracking of memory usage (RSS) and peak-resident sizes to validate `malloc_trim` benefits.

---

## Examples & useful commands 💡

- Run with pinned threads and CPU affinity (Linux):
  - `taskset -c 0-23 ./bin/cycle_finder --input graphs/huge_graph.bin --threads 24`

- Compare results between branches:
  - `git checkout main && make && mv bin/cycle_finder bin/cycle_finder.main`
  - `git checkout optimizations && make && mv bin/cycle_finder bin/cycle_finder.opt`
  - `./bin/cycle_finder.opt --input graphs/huge_graph.bin --threads 24 -o out.opt.json`
  - `./bin/cycle_finder.main --input graphs/huge_graph.bin --threads 1 -o out.main.json`
  - `python scripts/compare_results.py out.main.json out.opt.json`

---

## Changelog and rationale 📝

- 2025-12-29 — Introduced per-thread buffers + serial merging, lock-free visited bitmap, and per-thread pools to reduce contention and allocator pressure. Rationale: scale beyond ~24 cores and reduce per-edge overhead.

