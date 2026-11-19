# Array chimeras and the same-region evaluation

## What was wrong

A spacer's sequence occurs twice in the graph: in its CRISPR array and in the
phage it came from. At each end of the spacer the graph branches into the array
repeat or into the phage flank. Following the repeat (source-array crossing)
does not reach the target's flank. It walks along the array into the next spacer
and from there into that spacer's phage. The saved contig is phage DNA joined
through a CRISPR repeat to the target. It can pass the 80% identity / 50%
coverage criterion although the matched part is not around the target.

Evidence from the evidence19 outputs (32 paired COPSAC samples, local BLAST
against the Shah et al. catalogue; counts are FASTA records):

| Stage | Contigs | Contig matches | Target in the same reference region |
|---|---:|---:|---:|
| Discovery (seed) | 10,078 | 1,068 | 219 |
| Support 0.998 | 4,309 | 579 | 215 |

Among discovery contigs that crossed the source array, 41 of 675 matching
records had the target in the same reference region; in 525 of the 610 that
missed it, the sample's repeat lay between the matched sequence and the target.
The same pattern occurs without crossing when a flank enters another array that
uses a reported repeat.

## Changes (protospacer-workflow-v10, protospacer-seeded-local-v7)

- `--source-array-crossing` now defaults to `false` in both stages.
- New `--crispr-repeat-barrier true` (default): a flank stops at the first
  18-base word of any CRISPR repeat reported in the sample, consensus or
  per-spacer variant, either strand, including a word straddling the target
  boundary. The flank keeps the sequence before the repeat; stop reason
  `crispr_repeat`. An explicitly enabled crossing pass keeps its own rules.
- `run_summary.json` / `summary.json` record `crispr_repeat_barrier`,
  `crispr_repeat_word_length` and `crispr_repeat_words`.
- Old behaviour: `--source-array-crossing true --crispr-repeat-barrier false`.
- `scripts/copsac_batch.py` passes the new defaults and expects v10.

## Tests

- `mcaat phage-curate --self-test`: barrier on both strands and traversal
  directions, at the target boundary, not behind the target, disabled, and
  exempt for an explicit crossing pass.
- `tests/check_array_chimera.py BINARY OUT`: synthetic graph with a crossing
  chimera and a second-array chimera, both stages at support 0.998. With the
  defaults no saved extension contains a reported repeat and the reported
  neighbour still extends inside its own phage; without the barrier the
  second-array chimera returns. `--expect-chimera` confirms that an older binary
  produces both chimeras.
- `tests/check_seeded_extension.py` runs its source-passage case with crossing
  enabled explicitly and adds a default run that must keep repeats out of
  every extension.

## Same-region evaluation

`scripts/evaluate_target_region.py` keeps the contig criterion and adds: a
qualifying HSP covers at least 90% of the target interval, or the target itself
aligns to the same record and strand within 50 bp of the position the contig
HSP predicts. Remaining matches are classified by where their target lands.

```bash
python3 scripts/evaluate_target_region.py --root /vol/d/data/copsac \
  --db /path/all_curated_db --out region_eval --threads 16 \
  --stage old_seed=evidence19_seeds --stage old_support=evidence19_support_0.998 \
  --stage v10_seed=v10_seeds --stage v10_support=v10_support_0.998
```

Only samples with every listed stage complete are scored. `--paper-blast` uses
the manuscript's contig search (blastn-short, word 9, E 1000).

## Rerunning the COPSAC stages

`scripts/run_evidence19_stages.sh` repeats both evidence19 stages per sample
into new `<LABEL>_seeds` and `<LABEL>_support_0.998` folders with the same
parameters (`MCAAT=... LABEL=v10`). `BARRIER=false LABEL=v10nobarrier` gives the
crossing-off-only ablation.
