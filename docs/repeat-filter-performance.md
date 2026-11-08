# Repeat-filter runtime fix

The initial graph-based filter (`467eff2`) rebuilt the complete cycle-owner
index and materialized all candidate pairs after every duplicate or overlap
merge. Even unrelated systems were reprocessed. Dense duplicates also produced
quadratically many pairs, most of which were discarded after the first removal.

The replacement keeps stable system indices, cached cycle/role metadata, and a
persistent owner index. It updates only affected systems and visits candidates
in the same deterministic order. An unchanged retained system can remove its
remaining duplicates without restarting the comparison pass. Source IDs remain
sorted without sorting the entire accumulated list on every append. Graph-label
validation avoids temporary substrings and copies of cached k-mer labels.

## Local timings

Recorded on 2026-11-08, Intel Core i5-8259U, 16 GiB RAM, macOS, Apple Clang
with `-std=c++17 -O3`. Runs were sequential, with one recorded run per cell.
All cases use seed 42 and identical mode. Complete outputs match for both sizes
run against the old implementation.

| Candidate groups | Repeat families | Input cycles | Before (`467eff2`) | After |
|---:|---:|---:|---:|---:|
| 200 | 25 | 6,400 | 18.36 s | 0.62 s |
| 400 | 50 | 12,800 | 70.83 s | 1.23 s |
| 1,600 | 200 | 51,200 | Not run | 5.01 s |
| 2,000 | 1 | 64,000 | Not run | 7.24 s |

The 400-group case is approximately 58 times faster. Doubling from 200 to 400
groups increased the old runtime about 3.9-fold and the new runtime about
2-fold. This describes these duplicate-heavy fixtures, not every overlap pattern.

## Correctness checks

- 106 C++ filter checks, also run under AddressSanitizer and UndefinedBehaviorSanitizer.
- 360 seeded before/after comparisons: identical, overlapping/reverse-strand,
  repeat-variant, and ambiguous-repeat fixtures at 60%, 85%, and 100% support.
  Complete serialized systems match, including candidate cycles, observed
  sequences, boundaries, support, source IDs, alternatives, and original counts.
- Saved cycles from `CP065748.1_excerpt` and `CP003198.1`, at both 60% and 85%:
  all five output files match byte for byte against `467eff2`.
- Real-SDBG downstream integration: 5,911 checks passed, including edge/path
  validity, sequence reconstruction, deterministic reruns, stale output cleanup,
  rejection of damaged sidecars, and invalid representative edge slots.
- A dense-duplicate regression requires 256 systems to need exactly 255 pair
  checks and 256 initial index builds. Separate checks cover superset replacement
  and a merge that makes a previously disjoint group adjacent.

These establish agreement with the previous filter on the tested inputs, not
biological accuracy against a labeled reference.

## Reproduce

From the repository root, compile both implementations with the same compiler
and flags. The benchmark driver is identical; only the included filter changes.

```bash
mkdir -p build/repeat-filter-before/crispr_array
git show 467eff2:include/crispr_array/graph_repeat_filter.h \
  > build/repeat-filter-before/crispr_array/graph_repeat_filter.h
c++ -std=c++17 -O3 -I build/repeat-filter-before \
  bench/benchmark_repeat_filter.cpp -o build/benchmark-repeat-before
c++ -std=c++17 -O3 -I include \
  bench/benchmark_repeat_filter.cpp -o build/benchmark-repeat-after

build/benchmark-repeat-before 25 8 build/repeat-before.txt 42 identical
build/benchmark-repeat-after 25 8 build/repeat-after.txt 42 identical
cmp build/repeat-before.txt build/repeat-after.txt

python3 tests/compare_repeat_filter.py \
  build/benchmark-repeat-before build/benchmark-repeat-after \
  build/repeat-differential --seeds 30
```

The first two benchmark arguments are the number of independent repeat families
and the number of candidate groups per family. Each family has 32 distinct
repeat/spacer cycles; the identical mode duplicates those cycles across its
groups with different rotations. Optional modes generate partial overlaps,
reverse strands, sequence variants, and competing repeat regions. Fixture
construction and serialization are outside the reported `filter.run` timing.

The filter remains serial and still performs inference for each candidate group.
These benchmarks isolate the duplicate-filter scaling problem; they are not a
prediction of runtime for the complete 57-genome read dataset. Graph building,
cycle enumeration, and file I/O are outside these timings.
