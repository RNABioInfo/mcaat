# Sequencing evidence for protospacer extension

Research checked on 2026-11-16, before the experiments described in the
[follow-up evaluation](protospacer-evaluation-SRR14017553-v4.md). The proposed
coverage-ranking replacements failed to improve the benchmark and were not retained.
Only the smaller placement-continuation fix remains in the experimental workflow.
Prediction remains native to MCAAT and independent of the reference database.

## What was actually sequenced

The ENA experiment associated with **SRR14017553** is **SRX10394259**, in
**PRJNA715601**. Its submitted description specifies bulk faecal DNA extraction,
KAPA HyperPrep library preparation, and paired-end approximately 150-bp sequencing
on Illumina NovaSeq 6000. This is a bulk metagenome. The record does not document
the short MDA protocol used for the separately prepared viromes.
[ENA experiment record](https://www.ebi.ac.uk/ena/browser/view/SRX10394259).

Shah et al. (2023) constructed the 10,021-vOTU reference catalogue from enriched
viromes across 647 infants. Redgwell et al. (2025) describe the paired metagenome
and virome protocols explicitly: the bulk metagenomes used KAPA/NovaSeq, whereas
the viromes underwent a short multiple-displacement amplification step. Therefore
we must not attribute the graph's count variation to virome MDA without evidence.
The catalogue is also not a complete list of genomes present in this one sample.
[Shah et al.](https://www.nature.com/articles/s41564-023-01345-7),
[Redgwell et al., Methods](https://www.nature.com/articles/s41522-025-00674-1).

## What the studies establish

| Study | Relevant finding | Consequence for MCAAT |
| --- | --- | --- |
| [Roux et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5610896/) | In simulated viromes, populations at 2–5× depth often assembled well, but low depth and closely related strain mixtures impaired reconstruction. | Low multiplicity alone is not proof of error. Their genome-depth results are not a universal minimum edge-count threshold. |
| [Rihtman et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4893331/) | Approximately 30× depth sufficed for most isolate genomes; uneven experimental coverage motivated recommending higher depth for reliable isolate assembly. | Average depth and local gaps are different. These isolate results cannot be applied directly to a mixed, low-abundance metagenome. |
| [Sutton et al., 2019](https://link.springer.com/article/10.1186/s40168-019-0626-5) | Across simulated, mock, spiked, and gut viromes, low coverage and repeats impaired both contiguity and accuracy. | A longer graph path is not necessarily a better reconstruction. |
| [Haploflow, Fritz et al., 2021](https://link.springer.com/article/10.1186/s13059-021-02426-8) | Uses local, relative coverage to accommodate heterogeneous depth; stops when branch coverages cannot resolve the route. | Coverage is useful, but a global minimum or a single worst jump is not its complete method. Equal-support alternatives can remain unresolved. |
| [PhageTerm, Garneau et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557969/) | Collapsed direct terminal repeats can have roughly twice the sequence coverage of the rest of a phage. | A sharp increase can reflect a real repeated sequence. Counts alone do not identify a terminus; PhageTerm also uses read-start information. |
| [Phables, Mallawaarachchi et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10563150/) and [COBRA, Chen and Banfield, 2024](https://www.nature.com/articles/s41564-023-01598-2) | Combine graph structure and coverage with read-pair evidence to resolve or validate connections. | Their results do not establish that multiplicity alone can resolve every branch in our existing graph. We can study their principles without adding these tools to MCAAT. |
| [PenguiN, 2024](https://link.springer.com/article/10.1186/s40168-024-01904-y) | Read overlaps preserve linkage between variants that ordinary short-k de Bruijn graphs can lose. | Increasing beam width explores more routes but does not restore lost read linkage. |

Redgwell et al. additionally compared virome read recruitment inside predicted
prophages with their bacterial flanks in this same cohort. Their induction
analysis used paired datasets and a statistical test. It does not justify calling
a coverage jump in our bulk graph an induced prophage by itself.
[Redgwell et al., induction analysis](https://www.nature.com/articles/s41522-025-00674-1).

## Interpreting multiplicity correctly

The native workflow reads `EdgeMultiplicity` in `neighbors`. These are counts of
graph sequences, not a count of independent phages and not directly per-base
sequencing depth. Repeated or shared sequences pool contributions from every
origin that contains them. Matching spacer sequence can occur in the source
array as well as in a target genome, so its multiplicity is not automatically the
target phage's baseline.

A useful conceptual model for a sequence word w is

```
expected_count(w) = local_sampling_bias(w)
                  * sum_over_origins(abundance(origin) * copies(origin, w))
                  + contributions_from_sequencing_errors(w)
```

This is a proposed explanation of the observations, not an identifiable fitted
model: the graph does not label the origins, copy numbers, or error contributions.
For example, two genomes at depths 10 and 20 can contribute approximately 30 to
a shared segment. A route from a private segment of the first genome, through
the shared segment, and back to a private segment can legitimately look like
10 → 30 → 10. Compare its entry and exit; do not replace the baseline with 30.
This arithmetic is illustrative, not a measurement of SRR14017553.

Poisson sampling provides a starting null model for count variation, under
random fragment sampling. Real sequencing counts can be overdispersed; a
negative-binomial model permits variance larger than its mean. GenomeScope uses
this principle for k-mer spectra, but its diploid mixture model is not a model
of this viral community and should not be transplanted wholesale.
[Evans, Hower and Pachter, 2010](https://link.springer.com/article/10.1186/1471-2105-11-430),
[Vurture et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870704/).

Adjacent counts are dependent because the same reads supply overlapping words.
Likewise, one substitution can produce a run of erroneous k-mers rather than a
single isolated low count. Quake demonstrates this directly. Thus neither a
product of per-node probabilities nor a flat run of count-one nodes is calibrated
evidence for an independently supported long path. A low-count run can also be a
real rare variant, so its competing route, sequence difference and topology matter.
[Quake, Kelley et al., 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156955/),
[repeat-aware error modelling, 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3044310/).

## Problems identified in v3

The existing v3 implementation keeps the minimum unitig support for the entire
path and the maximum adjacent-window log-ratio. Both remember the worst past
event indefinitely. Its ranking prioritizes a band around the smallest such
maximum. Consequently, a locally smooth weak route can outrank a route with much
better downstream support after a short difficult region.

In the saved audit, the first reference-matching 750-base witness for Q000101 has
minimum count 1 and median count 10. Q000138 has minimum 1 and median 15. These
values are from `six_failure_traces.json`; they are reference-matching candidate
walks, not proof that every base or biological association is correct. Q000377
also shows why blindly favoring higher medians is insufficient: its matching
witness has median 1.

There is a separate confirmed selection problem: `process_query` stops testing
lower alignment tiers once any target is marked detected. Q000071, Q000101,
Q000138 and Q000229 have audited matching walks at placements that this rule
does not reach. Coverage rescoring alone cannot repair a placement never explored.

## Follow-up tests proposed from this evidence

1. Retain qualified alternative spacer placements within explicit work limits;
   report which placements were actually explored. A short exact placement must
   not prevent evaluation of a slightly mismatched placement with a supported
   extension.
2. Estimate local support and its variability on non-shared sequence. Inspect
   the distribution of counts, low-support run lengths, and transitions; stop
   using one lifetime minimum and maximum as the primary ranking criteria.
3. Explore through a shared region with a bounded beam and compare the private
   sequence before and after it. Keep competing hypotheses when both fit. High
   degree by itself does not prove a repeat.
4. Detect short diverge-and-rejoin alternatives and compare their spelled
   sequences and relative support. A low-support single-substitution bubble is
   different evidence from a long, distinct low-abundance branch. Neither is
   automatically an error.
5. Fit or calibrate any uncertainty model using graph observations and background
   controls, accounting for count dependence. Poisson or negative-binomial
   notation alone does not turn an unvalidated score into a probability of being
   correct. Do not substitute a fixed count cutoff for this calibration.
6. Keep adapter/source-array checks, readable rejection reasons, all original
   input observations and one primary prediction per query. Where evidence is
   insufficient, retain an explicit unresolved result rather than force a join.
7. Freeze the method, run all 398 distinct queries from the 530 observations and
   the same controls, then evaluate against the reference. Record gains and
   losses relative to both the published baseline and the rejected v3 run.

References can evaluate predictions and reveal defects during development; they
must not choose a path during prediction. This dataset is consequently a
development benchmark, not an independent estimate of generalization. The other
datasets should remain the subsequent external check.

The pre-experiment v3 baseline had 6 matching unique extensions out of 130 saved
at the agreed ≥90% identity / ≥50% query-coverage criterion. The final limited fix
saved 146 extensions with the same six matching regions, retaining their complete
sequences. One more input query supports an already detected region. This does not
establish an improved phage-extension method; see the full evaluation for lengths,
controls, failures and every matching sequence.
