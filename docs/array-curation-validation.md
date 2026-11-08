# Array curation repair: validation scope

The comparison baseline is `6205933`, before final candidate curation. The repair
adds configurable spacer-family/sequence-quality checks, exclusive ownership of
accepted circular sequences, complete review outputs for final curation, and an
immediate return when cycle enumeration reaches its cycle limit.

## Local saved graphs and annotated references

`bench/validate_array_references.py OLD_BINARY NEW_BINARY DATA_ROOT OUTPUT`
replayed all 24 available `*/cycles/str_paths.txt` fixtures against their actual
SDBG graphs. It read reference FASTAs from each fixture's `reference` directory,
canonicalized spacer sequences across both strands, and deduplicated identical
reference spacer sets **within each fixture**. Fixtures include excerpts and
repeated datasets; these totals are not independent genomes or a new biological
benchmark.

| Measurement, summed across fixture runs | Before | After |
|---|---:|---:|
| Reported candidates | 117 | 79 |
| Assigned cycle observations | 1,676 | 1,518 |
| Exact reference-spacer sequence matches | 1,031 | 1,031 |
| Reference sets with at least two exact spacers in one candidate | 66 | 66 |
| Circular sequences assigned to multiple candidates | 33 | 0 |

The reference input contains 133 deduplicated sets and 2,044 distinct spacer
sequences summed per fixture. **No previously matched reference spacer was
lost**, including when a non-representative variant provided the match. This
checks regression in sequence recovery. It does not establish locus precision,
complete biological recall, correct ordering, or all inferred boundaries.

The validator also checks one accepted owner per circular sequence, every
accepted spacer variant retained, representatives drawn from observed sequences,
and accounting of all final-curation observations in accepted or review files.

## Supplied metagenomic reports

The supplied reports contain 906 retained candidates and 24,009 assigned cycle
observations, with original support denominators. Replaying the **production
final-curation pass**, without reconstructing the original graph or rerunning
repeat inference, gives:

| Measurement | Result |
|---|---:|
| Accepted candidates | 319 |
| Accepted spacer families | 7,801 |
| Distinct accepted observed variants / cycle observations | 7,927 / 7,927 |
| Observations preserved in review | 16,082 |
| Missing or invented original observations | 0 |
| Circular sequences repeated in accepted output | 0 |

Curation of the 906 candidates took about 1.1 seconds locally, excluding report
loading/writing. This is not a timing of the complete server dataset. The
original 202,305-cycle input and its rejected groups were not available locally.
The replay preserves original unresolved-cycle counts in support denominators,
but cannot output unknown unresolved paths.

The concrete conflicting interpretations from original candidates 1, 474 and
537 fail the spacer-family similarity check. Original candidate 129 has 545
candidate cycles after merging and fails the exclusive 500-cycle limit. These
are explicit heuristic decisions, not ground-truth false-positive labels.

Build and reproduce a retained-candidate replay:

```bash
cmake -S tests -B build/check -DBUILD_SDBG_BENCHMARK=OFF
cmake --build build/check --target replay_array_candidates -j 2
build/check/replay_array_candidates /data/original-reports /data/curation-preview
python3 bench/check_array_curation_output.py /data/original-reports /data/curation-preview
```

`REPLAY_SCOPE.txt` records exactly what the helper did and maps candidate IDs to
original report IDs. Use the normal `postprocess-arrays` command for the complete
saved-cycle input, rather than treating this helper as full pipeline execution.

## Automated checks

- 2,507 candidate-quality checks, including 2,500 random comparisons against an
  independent dynamic-programming edit distance, complete-link variant grouping,
  observed variant preservation, minimum spacer families and tandem-like cases.
- 126 graph repeat-filter checks, including repeat-support thresholds, variants,
  boundaries, compatible merging, circular/strand ownership, residual support,
  single-residual-spacer rejection and preservation of competing assignments.
- Both suites pass AddressSanitizer and UndefinedBehaviorSanitizer.
- 37 CLI/pipeline cases: invalid parameters, settings, CLI overrides, and exact
  equality between an in-memory pipeline run and saved-cycle replay, including
  all new output files.
- Real-SDBG output integration checks validate actual downstream edge handles,
  sequence/path consistency, deterministic output and malformed-input rejection.
- `python3 tests/check_cycle_limit.py` extracts the production enumeration body
  into a tiny graph harness. At a limit of three, a three-cycle tangle exits in
  seven iterations and marks no nodes; with a limit of four all three cycles
  remain. A watchdog makes the former closing-edge loop fail promptly.

The default sequence thresholds are inspectable heuristics. All observed
variants and rejected final interpretations remain available for assessment;
changes in these counts must not be described as a measured improvement in
biological precision without full ground truth.
