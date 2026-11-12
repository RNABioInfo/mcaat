# SRR14017553: completed 81-input recovery check

**12 inputs have a connected, reference-consistent extension under the test below;
69 do not. No input remains unresolved under that test. Production selects a
matching extension for 11 inputs, with all ten previous matches preserved.**
The remaining recoverable input still loses at target-location selection.

The earlier claim that almost all 81 extensions could be recovered was unsupported.
The earlier 11-input diagnostic also omitted admissible target variants. This report
replaces that incomplete result with an explicit placement and continuation test.

[all81_results.tsv](all81_results.tsv) contains every original input, its yes/no
recovery result, the evidence, and the final production result.

## What was checked

The existing k=23 graph and all 81 frozen protospacers were reused unchanged.
The reference catalogue contains 1,041 previously identified locations in the
local Shah/EarlyVir database. For each input:

1. Find every connected graph target satisfying **90% substitution identity over
   at least 85% of the query**. There is no diagnostic placement/work cap. Short
   queries are checked at every offset inside a real graph k-mer. This produced
   **15,044 endpoint/context combinations across both directions**. These are
   graph search states, not counts of independent organisms or target loci.
2. From those endpoints, test every known reference-flank length from **100 to
   750 bases**, anchored immediately beside the target, allowing substitutions
   and indels at **at least 90% aligned identity**. Every emitted base must follow
   an actual graph edge. Reference ends and unknown bases end the known interval.
3. Exclude default adapters and closed walks entirely inside the source-array
   cycles. External revisits are permitted in the existence check. Production
   multiplicity gates and beam pruning are absent from this diagnostic.
4. Accept a positive only with a graph-validated path and independently checked
   alignment. Accept a negative only after exhausting all admissible starts and
   all eligible reference sides for that input. Interrupted/capped searches are
   never negative evidence; difficult cases received further completed searches.

The result is limited to this reference-location catalogue and the stated
ungapped **target-placement** rule. Target indels are outside that rule; indels
are allowed in the **flank** comparison. A negative does not label the original
protospacer false or assert anything about unknown reference sequences.

The evidence is saved in [audit_summary.json](audit_summary.json),
[complete_placements.tsv.gz](complete_placements.tsv.gz),
[connected_path_witnesses.tsv.gz](connected_path_witnesses.tsv.gz), and
[completed_negative_searches.tsv.gz](completed_negative_searches.tsv.gz).
`bench/collect_complete_recovery_audit.py` checks complete endpoint coverage,
alignment operations, and exhausted-search flags; its normal mode fails if an
input remains unresolved.

## Production comparison

All runs used 24 threads, the existing graph, a maximum of 750 bases per flank,
and at most three visits to a node outside source-array cycles. Reference DNA
was used only for separate diagnosis and evaluation, never production selection.

| Run | Inputs retained | Inputs extended | Matching extensions | Previous matches lost |
|---|---:|---:|---:|---:|
| Previous N=32 | 81 | 56 | 10 | — |
| New N=16 | 81 | 57 | 11 | 0 |
| New N=32 | 81 | 59 | 11 | 0 |
| New N=64 | 81 | 59 | 11 | 0 |

Here “matching” means local BLAST at **≥90% identity and ≥50% single-HSP query
coverage of the complete emitted extension**, as requested. This is a separate
measure from the anchored 100-base existence test. All previously emitted
predictions—not just their reference matches—remain sequence-identical at each
corresponding beam width. See [beam_comparison.tsv](beam_comparison.tsv) and
[previous_prediction_preservation.json](previous_prediction_preservation.json).
Other emitted extensions remain unvalidated predictions.

## The newly repaired failure

`seed_edac35e71117a388` has input `AGGAGGGGGCTCTGCCGCACGGCGC`. The observed graph
target `AGCAGGGGGCTCTGCTGCACCGCGC` satisfies the accepted 90%/85% interval rule,
but every full 23-mer anchor has more than one substitution. The previous fast
mapper therefore missed it. A bounded masked BOSS lookup now finds it after the
earlier placement attempts fail, sharing their existing work budgets.

With the target found, the original score selected a 42-base left flank. It took
the competing branch eight bases from the target; the reference-consistent branch
has count one there. The new rescue score tolerates ordinary low-depth count
variation and selects the verified **136-base left flank** instead. The complete
prediction is **310 bases**, including a 149-base right flank and the 25-base
observed target. It matches **OTU_6330 at 90.820% identity over 98.387% of the
prediction**. The original input FASTA record remains unchanged.

Applying this scoring change globally lost one previous match. That experiment
is recorded in [rejected_global_noise_comparison.tsv](rejected_global_noise_comparison.tsv).
The shipped policy applies it only to the final masked-target rescue, preserving
existing exact and first-fallback predictions. Source-array, adapter, coverage,
visit, and work limits still apply. The score is a heuristic, not a probability
of viral origin.

## One selection failure remains

`seed_12839d3c8ccafc84` is the 17-base input `ACAATAATTATGTTTTT`. A matching path
exists: isolating its observed target yields a **1,070-base exact match to
OTU_12626 at N=32**. Ordinary production still selects another admissible target
with a better multiplicity score. The matching candidate is present in the
placement diagnostics; it is not silently reported as absent.

The complete search found 5,132 endpoint/context combinations on each side of
this input, including overlapping k-mer contexts. The selected and matching
extensions both have unique 23-mer paths, so rejecting repeated-node walks would
not solve this choice. The additional source-spacer sequence agrees with only
5 of 13 adjacent bases of the matching path, so it supplies no strong extra
anchor either. **This target-location choice is not fixed.** Selecting it using
the database would violate the requested reference-independent production method.
The earlier [isolated-path and multiplicity evidence](../protospacer-target-placement-2026-11-12/README.md#the-one-remaining-selection-failure)
remains applicable.

## Validation and rerun

Release and ASan/UBSan tests pass. Tests compare the complete mapper with
exhaustive small-graph searches; diagnostic alignments with independent string DP;
and native BOSS lookup with full scans of small graph fixtures. Native CLI tests
verify the masked-target case, unchanged inputs, failure handling, and identical
1-/24-thread predictions. Full-graph N=16/32/64 runs verify all 81 sequences,
coordinates, and unchanged graph files. Checksums and logs are in
[validation.json](validation.json).

Run [rerun.sh](rerun.sh) on the server for the full pull/build/extension command
using the previously supplied paths and 24 threads. It writes to a new output
directory and reuses the graph. The original 81 inputs are included in the repo's
`bench/fixtures/protospacer-baseline-81/Protospacers.fasta`.
