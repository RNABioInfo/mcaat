# Remove redundant circular sequence sets

Conservative curation previously kept complete candidate subsets when graph-role
checks had prevented an earlier merge. Changing the repeat boundary changes the
reported spacer, so comparing repeat strings or spacer strings alone misses this.

The final curation pass now works on exact complete strings:

```
key(unit) = lexicographically smallest rotation of unit or reverse_complement(unit)
evidence(candidate) = set(key(observed_repeat + observed_spacer))
remove A if evidence(A) is a subset of evidence(B) for ONE retained B
```

Candidates are processed by decreasing evidence-set size. Equal-size sets are
ordered by higher repeat-path support fraction, longer eligible repeat, canonical
repeat sequence, then original candidate ID. The selected repeat and spacers
remain observed strings. The tie break is deterministic; it does not establish
that one boundary is the true genomic boundary.

An inverted index of already retained candidates is searched using the candidate's
rarest indexed circular sequence. Full subset checks are performed only against
that posting list. No global pairwise overlap matrix is allocated. If any sequence
has no retained owner, the candidate necessarily adds evidence and is retained.
Coverage by the union of multiple candidates is insufficient for removal.

This operation preserves the union of accepted circular sequences. It does not
preserve the number of displayed spacer strings: two cuts of the same circle can
give different strings. Source start IDs are merged into the surviving candidate,
without increasing its support count. Removed candidates receive the decision
`redundant_circular_subset`, linked to the surviving candidate ID. Original
repeat/spacer assignments remain in `CRISPR_graph_alternatives.tsv`; unresolved
source cycles remain in the review table when those paths were available. A source
cycle without an assigned valid repeat/spacer is not counted as an array observation
for the subset test. The candidate's original counts are retained in the decision.

## Length checks

The defaults are repeat length 23–55 nt and spacer length 23–50 nt. These inclusive
limits were already applied during repeat placement. Final curation now checks
them too, including report replay. Invalid assignments are excluded and retained
as alternatives; their cycles are not silently removed from the support denominator.
Support and the minimum two observed spacers are checked after exclusion. The four
configured bounds and exclusion counts are included in the output summary.

## Supplied 47–100 run: deduplication-only comparison

These measurements isolate the exact-set deduplication change. Current MCAAT also
runs native string filtering by default; see [string filtering](string-structure-audit.md)
for the combined results. Use `--string-filter off` only to reproduce this older
deduplication-only comparison.

| Measurement | Before | After |
|---|---:|---:|
| Accepted candidate reports | 934 | 904 |
| Reported spacer entries | 23,885 | 22,918 |
| Distinct complete circular sequences | 17,247 | 17,247 |
| Distinct spacer strings, merging reverse complements | 20,922 | 20,493 |
| Exact spacer matches to the available local references | 842 | 842 |

All reported repeat lengths were 23–55 and all spacers 23–50; none was outside the
default bounds. The new sequence-set pass removed 30 redundant reports and 967
accepted observation entries. All 23,885 original observations are accounted for
exactly across accepted and alternative tables. There are no remaining accepted
set-containment pairs. This is exact redundancy removal, independent of the
separate exploratory 74-candidate motif/conservation screen.

Final curation of the supplied reports took 2.56 seconds and used 216 indexed
subset comparisons. This timing excludes input parsing, output writing, graph
loading, cycle search and repeat inference. No MAFFT alignment is used here.

For example, all 219 assigned circles of original `Array_349` occur in original
`Array_127`, even though their different repeat boundaries give no exact spacer
matches. The smaller report is removed; its boundary alternatives remain available.
Array IDs are renumbered after filtering, so use `candidate_id` or the replay map
to trace original names.

## Reproduce on existing reports

```bash
cmake -S tests -B build/array-tools -DBUILD_SDBG_BENCHMARK=OFF
cmake --build build/array-tools --target replay_array_candidates -j 2
build/array-tools/replay_array_candidates /path/to/reports /path/to/new-output \
  --array-filter-mode conservative --string-filter off --cycle-min-length 47 --cycle-max-length 100 \
  --min-repeat-len 23 --max-repeat-len 55 --min-spacer-len 23 --max-spacer-len 50
python3 bench/check_array_curation_output.py /path/to/reports /path/to/new-output
```

The replay tool checks consistency of supplied paths and strings, recovers graph
labels from those observations, and calls the same final curation/output code as
MCAAT. It does not rerun inference or recover earlier rejected cycles. The source
report's original unresolved-cycle counts remain in its support denominators;
unavailable unresolved paths cannot be reconstructed. `REPLAY_SCOPE.txt` records
that limitation and the original array-ID mapping. Full `postprocess-arrays`
execution on the graph and saved cycles also applies the new pass automatically.

The regression suite covers shifted boundaries, reverse-complement rotations,
equal-set tie breaking, idempotence, partial overlaps, shared repeats with distinct
cycles, preservation of unresolved source paths, and a 256-report duplicate case
requiring only 255 subset comparisons. Boundary tests cover inclusive length
limits, invalid spacer assignments, and support-denominator preservation.

Validation passed 153 graph-filter checks (also under ASan/UBSan), 2,509 quality
checks, and 60 CLI/pipeline cases. Across 24 local saved-cycle reference fixtures,
all 1,105 previously matched reference spacers remained (summed across fixtures,
which include repeated genomes/excerpts). A second pass on the cleaned supplied
reports removed nothing; accepted paths, variants and text reports were byte
identical. Complete approximately 750-array locus-level recall was not measured.
