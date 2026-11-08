# String evidence after cycle detection

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

## Method

1. Keep each observed unit as `repeat | spacer`. Align the repeat and spacer
   regions separately with MAFFT (`--retree 1 --maxiterate 0`). Check that removing
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

## Separate review flags

- **Short-motif repetition:** at least two existing similarity families, with
  at least 85% of family-balanced observation weight fitting a motif at 90%
  identity or higher.
- **Conserved spacer:** at least three families, at least 80% of spacer alignment
  columns conserved at 85% or higher, and mean repeat minus spacer conservation
  no greater than 0.15.
- **Boundary ambiguity:** at least three consecutive conserved columns at a
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

## Supplied 47–100 run, 2026-11-08

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

These are reasons to inspect strings, not newly established biological labels.
Neither rule is installed as a production rejection default. Full annotations
and separate validation data are needed before choosing a production policy.

## Tests and methodological context

Run `python3 tests/test_string_audit.py`. Twelve tests cover family balancing,
gap handling, motif substitutions, reference-set deduplication, invalid inputs,
preservation of strings by real MAFFT, boundary flags, and an end-to-end fixture
where tandem-like strings are flagged while all positive-control spacers remain.
MAFFT-dependent tests explicitly skip if MAFFT is absent.

[CRISPRDetect](https://doi.org/10.1186/s12864-016-2627-0) uses alignments to refine
repeat boundaries and distinguish tandem repeats. This audit uses that general
principle, not its exact ClustalW procedure or published score.
[CRISPRidentify](https://doi.org/10.1093/nar/gkaa1158) reports overlap in spacer
similarity between positive and negative examples. That supports keeping a
single similarity threshold from acting as sufficient biological evidence.
