# String evidence after cycle detection

## Automatic filtering in MCAAT

Both **ordinary `mcaat` execution and `mcaat postprocess-arrays` now apply the
string filters by default**. No separate audit, Python, MAFFT installation or
extra command is required. The implementation uses MCAAT's bundled SPOA library
for global multiple sequence alignment of observed strings; it never generates
a consensus or changes a repeat/spacer boundary. The resulting filtered systems
feed the normal array reports and downstream CAS/phage stages.

Filtering follows structural and length checks and precedes final circular-set
deduplication. It rejects candidates meeting either the short-motif or
conserved-spacer rule below. Every excluded observation remains in
`CRISPR_graph_alternatives.tsv`; `CRISPR_graph_decisions.tsv` records the reason.
`CRISPR_string_scores.tsv` records the numerical evidence for every candidate,
with `NA` for scores not computed. The summary and normal run's `parameters.json`
record the active settings. These are computational heuristics, not calibrated
probabilities or proof that every remaining candidate is a CRISPR array.

| CLI option | Settings-file key | Default |
|---|---|---:|
| `--string-filter` | `string-filter` | `on` |
| `--string-column-support` | `string-column-support` | 0.85 |
| `--string-motif-identity` | `string-motif-identity` | 0.90 |
| `--string-motif-support` | `string-motif-support` | 0.85 |
| `--string-conserved-spacer-fraction` | `string-conserved-spacer-fraction` | 0.80 |
| `--string-max-contrast` | `string-max-contrast` | 0.15 |

`--string-filter off` is an explicit comparison override. Selecting conservative
array-filter mode does **not** disable these two new rules; it controls the older
pairwise-similarity and cross-candidate ownership heuristics. At the existing
cycle cap, native alignment is skipped to bound work and the candidate gets
`string_alignment_skipped_at_cycle_limit`; motif scoring still runs. Skipped
alignment is counted explicitly in the summary.

Native replay of the supplied 934-candidate reports rejects **71** candidates by
string evidence and **24** redundant reports, leaving **839 candidates / 16,776
reported spacer entries**. The 7,109 excluded observations remain in alternatives.
All **842** exact spacer matches to the available local reference subset remain.
The measured C++ curation time was 3–7 seconds locally (the latter during compilation);
it is not the runtime for graph loading or cycle search.

Native SPOA alignment and the earlier MAFFT audit do not give identical column
profiles. MAFFT additionally flags `Array_132`, `Array_354` and `Array_777`
(77 entries together); the native pass has no additional rejections relative to
that audit. Native alignment uses lexical ordering of unique observed strings,
global alignment scores `(match=3, mismatch=-5, gap=-3)`, and then restores all
duplicate observations for family-balanced scoring. Neither method substitutes
the aligned strings for the observed output sequences.

## Optional standalone MAFFT audit

`bench/audit_string_structure.py` audits the observed repeat/spacer strings in an
existing MCAAT output directory. It needs Python 3 and MAFFT. It does not load the
graph, enumerate cycles, change boundaries, reconstruct a consensus sequence,
or remove production candidates. Its numerical cutoffs are exploratory review
settings, not a calibrated CRISPR classifier.

```bash
python3 bench/audit_string_structure.py /path/to/mcaat-output results/string-audit \
  --workers 2
```

The input directory must contain `CRISPR_graph_summary.tsv`,
`CRISPR_graph_systems.tsv`, `CRISPR_graph_paths.tsv` and
`CRISPR_spacer_variants.tsv`. Use a fresh output directory. Each worker invokes
MAFFT with one thread. The optional `--references /path/to/data` reads spacer
records from `data/*/reference/*`; FASTA headers must start with `spacer`.

### Reproduce the earlier MAFFT-based filtering

The optional audit itself only scores candidates. To reproduce its exact
MAFFT-based decisions, rather than native SPOA decisions, use the companion step:

```bash
cmake -S tests -B build/array-tools -DBUILD_SDBG_BENCHMARK=OFF
cmake --build build/array-tools --target replay_array_candidates -j 2
python3 bench/apply_string_filter.py /path/to/mcaat-output results/string-audit \
  results/string-filtered --replay-binary build/array-tools/replay_array_candidates
```

Use the exact input directory that was audited, and a new output directory. The
application verifies SHA-256 hashes of all four input reports before using the
audit, checks complete candidate-ID coverage, and rejects stale audits. It applies
the short-motif/conserved-spacer union before exact circular-set deduplication.
Boundary-only and single-family flags do not cause these string rejections.
The companion explicitly disables the native pass to avoid mixing two aligners'
decisions; it is not needed for production MCAAT execution.
The supplied original cycle range and repeat/spacer length bounds are retained
(older reports lacking repeat/spacer metadata use the existing 23–55/23–50
defaults). No graph, cycle search or new alignment is needed for an existing audit.

All removed observations are written to `CRISPR_graph_alternatives.tsv` with
`string_filter:<reason>` or `redundant_circular_subset` decisions. The manifest
`string_filter_rejections.tsv` uses original array IDs. Accepted reports receive
new array IDs; `candidate_id` and `REPLAY_SCOPE.txt` preserve the mapping.
`string_filter_summary.json` records parameters, input/audit hashes, counts, and
available exact-reference comparison results.

On the supplied 934-candidate run, applying the existing thresholds rejects 74
candidates containing 6,609 spacer entries. Deduplicating the remaining evidence
removes 24 more reports containing 577 entries. The result is **836 candidates
and 16,699 reported spacer entries**, down from 934 and 23,885. All 7,186 excluded
observations remain in alternatives. All 842 previously detected exact matches
to the available local reference spacers remain. This is not a measurement of
precision/recall against the full approximately 750 annotated arrays.

The historical audit thresholds remain configurable in the audit command.
They were not adjusted to obtain a requested number of candidates.

## Method

1. Keep each observed unit as `repeat | spacer`. Align the repeat and spacer
   regions separately with native SPOA (or MAFFT `--retree 1 --maxiterate 0`
   in the standalone audit). Check that removing
   alignment gaps recovers every original string exactly.
2. Give each existing spacer similarity family equal total weight. Divide its
   weight equally among its observations. This limits domination by a family
   with many variants; these families are not asserted to be biological loci or
   independent spacers.
3. At each alignment column, calculate the largest weighted A/C/G/T fraction.
   Gaps stay in the denominator and cannot supply conserved-base evidence.
4. Test each unaligned repeat+spacer unit for short-motif repetition. For every
   period from 1 through floor(unit length / 3), fit the most frequent base at
   each position modulo that period. Record the best matching-base fraction.
   This permits substitutions, requires at least three nominal motif copies,
   and does not correct indels. Different units can select different motifs.

The motif fit is a descriptive identity, not a statistical significance or
CRISPR probability. Conservation is alignment-dependent. The supplied graph
repeat is identical across observations in every candidate in this run, so its
100% conservation is imposed by extraction, not independent biological evidence.

## String rules and separate review flags

- **Short-motif repetition:** at least two existing similarity families, with
  at least 85% of family-balanced observation weight fitting a motif at 90%
  identity or higher.
- **Conserved spacer:** at least three families, at least 80% of spacer alignment
  columns conserved at 85% or higher, and mean repeat minus spacer conservation
  no greater than 0.15.
- **Boundary ambiguity (standalone audit only):** at least three consecutive conserved columns at a
  spacer end, leaving at least five interior columns with mean conservation
  below 85%. This is not included in the preceding two flags' union and is not
  an instruction to discard the candidate or automatically extend its repeat.
- **Single similarity family:** reported separately; not an automatic rejection.

All cutoffs are configurable; see `--help`. The audit's column and motif support
settings are separate from MCAAT's graph `--repeat-support` parameter.

`candidate_scores.tsv` contains the scores and flags. `column_profiles.tsv` and
`alignments/Array_N.txt` expose the positional evidence. `rule_impact.tsv` shows
the hypothetical effect of excluding flagged candidates on exact local-reference
spacer matches. Matching is reverse-complement invariant. Reference spacer sets
are deduplicated by sequence set, so reference counts are not genomic locus counts.
`summary.json` records parameters, MAFFT version, runtime and SHA-256 input hashes.
Inputs are hashed again after the audit to catch changes during processing.

## Historical MAFFT audit of the supplied 47–100 run, 2026-11-08

| Exploratory rule | Flagged candidates | Reported spacer entries in those candidates |
|---|---:|---:|
| Short-motif repetition | 55 | 6,352 |
| Conserved spacer | 21 | 305 |
| Union (two candidates overlap) | 74 | 6,609 |

There are 934 input candidates and 23,885 reported spacer entries. Hypothetical
exclusion of the union would leave 860 candidates and 17,276 entries. Shared
spacers can occur in multiple candidates: the number of globally distinct
spacer sequences, counting reverse complements once, would change from 20,922
to 16,767. Nine candidates have the separate boundary flag.

The local positive controls contain 129 distinct reference spacer sets. Across
the supplied predictions, 842 distinct spacers match these controls exactly;
47 candidates contain at least two spacers from a single reference set. None of
the 74 flagged candidates contains an exact local-reference spacer match, and
all 842 detected matches would remain. This does not establish recall against
the full approximately 750 annotated arrays or precision of the remaining 860.
An unmatched candidate is not thereby a false positive.

Examples in this run:

- `Array_29`: one 63-base unit fits three copies of a 21-base motif with 61/63
  matching bases; the other fits with 59/63. This is direct periodicity evidence.
- `Array_43`: only 6 of 49 spacer alignment columns fall below 85% conservation.
- `Array_86`: 11 conserved columns at the right spacer end suggest examining the
  reported repeat boundary, while the spacer interior remains variable.
- `Array_356`: variable spacers and 85 exact spacer matches to one local reference;
  neither proposed rejection-style rule flags it.

These are computational rejection rules, not established biological labels.
The native versions now run by default as described above. Full annotations and
separate validation data are needed to measure their classification accuracy.

## Tests and methodological context

Build the standalone C++ tests and run `ctest --test-dir build/array-tools
--output-on-failure`. `test_string_structure` checks motif scores against
exhaustive templates, alignment preservation and order independence, family
weighting, gap handling, the two rules, bounded alignment and option validation.
The graph-filter suite proves default rejection from actual simple cycles and
restoration with the explicit off override. `tests/check_graph_pipeline.py`
checks the compiled production binary, normal versus saved-cycle output equality,
default scoring, settings and CLI overrides.

Validation of the integrated default passed 140 native string checks, 163 graph
filter checks, 2,509 candidate-quality checks, and 88 compiled CLI/pipeline cases.
The graph and string checks also passed ASan/UBSan. Across 24 saved-cycle reference
fixtures, retained candidates changed from 114 to 110 and observations from 1,757
to 1,741, while all 1,105 existing exact reference-spacer matches remained. These
matches are summed across fixtures, which include repeated genomes/excerpts.
A fresh cycle search through ordinary `mcaat --graph ...` on the local
`CP071793.1` fixture scored 61 candidates and automatically rejected five by the
native string rules, without a separate audit or filter invocation.

Run `python3 tests/test_string_audit.py`. Twelve tests cover family balancing,
gap handling, motif substitutions, reference-set deduplication, invalid inputs,
preservation of strings by real MAFFT, boundary flags, and an end-to-end fixture
where tandem-like strings are flagged while all positive-control spacers remain.
MAFFT-dependent tests explicitly skip if MAFFT is absent.

`python3 tests/test_apply_string_filter.py --replay-binary build/array-tools/replay_array_candidates`
checks stale-input rejection, candidate-ID coverage, flag consistency, exclusion
of boundary-only flags, CRLF inputs, and end-to-end preservation of rejected
observations and positive-control spacers. The C++ graph-filter suite also checks
that a string-rejected parent cannot suppress an unflagged contained candidate.

[CRISPRDetect](https://doi.org/10.1186/s12864-016-2627-0) uses alignments to refine
repeat boundaries and distinguish tandem repeats. This audit uses that general
principle, not its exact ClustalW procedure or published score.
[CRISPRidentify](https://doi.org/10.1093/nar/gkaa1158) reports overlap in spacer
similarity between positive and negative examples. That supports keeping a
single similarity threshold from acting as sufficient biological evidence.
