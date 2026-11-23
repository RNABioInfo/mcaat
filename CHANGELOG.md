# Changelog

## Protospacer boundary walk (tools/boundary_walk)
- New `tools/boundary_walk.cpp`: `FindProtospacers` and `ProtospacerExtend` from the thesis, graph topology and multiplicity only, with seven exchangeable beam-search criteria. Both boundary sets are still required to be non-empty.
- The nodes of a spacer are no longer required to form one connected path, which the thesis never asks for. On SRR14017565 requiring it discarded 621 of 1,962 spacer observations; dropping it took the two-sided candidates from 245 to 575 and the strict matches from 22 to 93.
- The multiplicity window is keyed to the running median of the walk instead of to the boundary node.
- A walk is kept when the harmonic mean of its node multiplicities, pooled over both flanks, is at least `--min-harmonic-mean` (default 3).
- On SRR14017565 that gives 268 contigs with 93 strict matches at 56 loci, against 286 contigs with 54 matches at 5 loci from `mcaat phage-curate`. 52 of the matches have the alignment covering the spacer and at least 50 bases on each side; the production run has none. One sample only.
- New `tools/neighbour_check.cpp`: checks a neighbour lookup against `SDBG::IncomingEdges`/`OutgoingEdges` on every k-mer of every spacer.
- New `scripts/run_boundary_walk.sh`, `scripts/eval_criteria.py`, `scripts/rule_table.py`, `scripts/spacer_check.py`. See `docs/boundary-walk.md`.
- `boundary_walk` walks the spacers in parallel, `--threads N` to limit it. Output is written in source order and is byte-identical to the serial version on SRR14017565; 16 cores took that sample from 29m36s to 8m20s, the remainder being the serial graph load.


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
