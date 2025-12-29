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

If you'd like, I can also:
- Add a short benchmark script that automates runs at multiple thread counts and produces a small CSV/plot; or
- Add micro-benchmarks for the visited bitmap vs a lock-based approach for small stress tests.

