# boundary_walk: the thesis extension, with the spacer no longer required to be one path

`tools/boundary_walk.cpp` implements `FindProtospacers` and `ProtospacerExtend`
as the thesis defines them, with the multiplicity criterion of the beam search
exchangeable, and reads nothing but graph topology and multiplicity. No
reference is used to build or to select anything.

For every array the tangle is the node set of its repeat and its spacers. For
every node of a spacer the neighbours outside the tangle are collected: those
pointing in give B_in, those pointed to give B_out. A spacer is a candidate only
when **both** sets are non-empty. Each node of B_out is extended forward and
each node of B_in backward, beam width 50, up to 1250 added bases, admitting a
node whose multiplicity lies inside a window and lies outside the tangle. The
best path of each side is kept and the record is `left + spacer + right`.

## What differs from the thesis

**The window is keyed to the running median of the walk, not to the boundary
node.** The thesis window is `[0.1*mult(b), 5*mult(b)]` for the boundary node
`b` the walk started from. A boundary node that happens to sit on a
low-coverage k-mer then sets a window the rest of the flank cannot stay inside.

**The nodes of a spacer are not required to form one connected path.** The
thesis asks for the neighbours of *every node v of a spacer*; it never asks
those nodes to form a path. Requiring one discarded 621 of 1,962 spacer
observations on SRR14017565. Of the candidates that requirement was throwing
away, 178 have every k-mer of the spacer present in the graph and only one
internal edge missing. Both boundary sets are still required to be non-empty,
so the two-sided definition is unchanged.

A missing internal edge does create boundary nodes in the middle of a spacer,
which could in principle satisfy the two-sided test without either side being a
real flank. Measured, the opposite happens: of the recovered matches 44 of 71
have the alignment reaching both flanks, against 8 of 22 for the candidates
whose spacer did form one path.

## The rule

The criterion chooses the path; it does not decide whether to keep it. A walk is
kept when the harmonic mean of its node multiplicities, pooled over both flanks,
is at least `--min-harmonic-mean`. A walk whose harmonic mean is 1 was built
entirely out of k-mers seen once, which a sequencing error also produces.

On SRR14017565 every strict match survives up to 3 and all but one up to 5:

| `--min-harmonic-mean` | contigs kept | strict matches | rate | loci |
|---:|---:|---:|---:|---:|
| 0 | 575 | 93 | 16.2% | 56 |
| 2 | 290 | 93 | 32.1% | 56 |
| **3** | **268** | **93** | **34.7%** | **56** |
| 5 | 242 | 92 | 38.0% | 56 |
| 8 | 185 | 66 | 35.7% | 44 |

3 is the default because it is inside that plateau rather than on its edge. The
threshold was read off this curve, on this one sample.

## SRR14017565

A *strict match* is one alignment with at least 80% identity over at least 400
aligned bases whose aligned part covers at least 90% of the spacer. The spacer's
position comes from the FASTA header, not from the hit, so a flank matching
something else does not count. Two alignments on one reference record are the
same locus only if their starts are within 100 bases. Local BLAST against the
Shah et al. catalogue.

The production run compared against is the 286-contig `mcaat phage-curate` run
on this sample, produced before the repeat barrier existed: its recorded config
has `source_array_crossing: true` and no barrier. Reproducing that baseline from
this branch needs `--source-array-crossing true --crispr-repeat-barrier false`,
because both defaults changed afterwards.

| Run | Contigs | Strict matches | Rate | Loci | Alignment spans both flanks |
|---|---:|---:|---:|---:|---:|
| `mcaat phage-curate` (286-contig run) | 286 | 54 | 18.9% | 5 | **0** |
| boundary_walk, harmonic_mean, rule at 3 | 268 | 93 | 34.7% | 56 | not counted |
| boundary_walk, harmonic_mean, rule at 5 | 242 | 92 | 38.0% | 56 | 52 |

Fewer contigs, more matches. The last column counts matches where the alignment
covers the spacer and at least 50 bases on each side of it, which is the direct
test of the two-sided definition; it was counted on the rule-at-5 run, so the
default's figure is at least 52 and was not measured separately. Every one of the production run's 286 contigs
is extended on one side only — 172 left, 114 right, 235 of them by the
pipeline's own choice of the alternative flank and 51 because the other side
could not be extended — so all 54 of its matches reach one flank and none
reaches two.

The seven criteria on the same candidate set, at the default rule. The rule is
applied after the criterion has chosen its path, so each criterion keeps a
different number of contigs:

| Criterion | Contigs kept | Strict matches | Rate | Loci |
|---|---:|---:|---:|---:|
| harmonic_mean | 268 | 93 | 34.7% | 56 |
| geometric_mean | 268 | 77 | 28.7% | 49 |
| mean_multiplicity (thesis) | 264 | 69 | 26.1% | 47 |
| length_weighted_mean | 264 | 69 | 26.1% | 47 |
| boundary_stable | 163 | 52 | 31.9% | 33 |
| bottleneck | 124 | 50 | 40.3% | 37 |
| locally_smooth | 91 | 29 | 31.9% | 25 |

`bottleneck` has the better rate — 40.3% against 34.7% — on a third of the
contigs, and would be the choice if the rate were the only thing being asked
for. It is not the choice here, because 50 matches is fewer than the production
run's 54: it gives up the "more matches" half of the goal. `harmonic_mean` is
the only criterion that is below the production run on contigs and above it on
matches at the same time.

Match rate rises with how far the walk got: 0% below 400 bp of flank, 31% at
the 1250 bp limit on one side, 60% when both sides reach it. Of the candidates,
165 never exceed 150 bp on either side and cannot match at all.

## Negative result: the window's lower bound

88% of all refused steps are refused for falling below the window, which looked
like the thing holding the walks short. Lowering the bound from 0.1 to 0.02,
0.002 and 0 makes the result worse, not better — on the connected-only set the
matches fall from 22 to 15 and the loci from 16 to 13. The bound is keeping the
walk off low-coverage error branches. It is left at 0.1.

## Not settled

- One sample. SRR14017570 is the sharper test, because nothing it produces
  matched in the earlier runs.
- The threshold was chosen from its own curve.
- 57 candidates have under half their k-mers in the graph and still match at
  26.3%, which is 15 of the 93 matches. There is no mechanism for that yet, so
  that slice is a reason for doubt rather than for confidence.
- boundary_walk and the production run share only 2 loci, out of 56 and 5.
- The tangle blocks the source array only. 8 of 389 contigs from the earlier
  production run crossed a *different* array's repeat.

## Running it

On the machine with the graph:

```bash
BW=/path/to/build/boundary_walk bash scripts/run_boundary_walk.sh SRR14017570
```

Every setting is left at the binary's defaults, which are the configuration
above. `EXTRA="--min-harmonic-mean 5"` and the other flags override them for an
ablation; `boundary_walk` with no arguments prints the defaults.

Scoring needs the reference database, not the graph:

```bash
python3 scripts/eval_criteria.py --dir SRR14017570/bw --db all_curated_db --out eval570
python3 scripts/rule_table.py --dir SRR14017570/bw --eval eval570 --criterion harmonic_mean
```

`eval_criteria.py` reports every criterion under several match definitions;
`rule_table.py` applies the recorded harmonic mean and scores what it keeps.
`spacer_check.py` reports how many of the long alignments actually cover the
spacer rather than a flank.
