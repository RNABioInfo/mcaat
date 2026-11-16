# Adapter-context repair: SRR14017553

The v6 discovery check could accept a technical adapter branch as evidence that a spacer had a location outside its source array. This could then authorize source-crossing extension along a different branch. Checking only the final extension FASTA would miss this error.

## What failed

The old filter recognized two longer TruSeq cores with at most one substitution. Twenty-one accepted targets had a selected discovery context containing the exact shared prefix `AGATCGGAAGAGC`, followed by sequence that escaped those longer-core tests. The 21 targets were associated with 20 unique extensions; neither those targets nor those extensions passed the fixed reference criteria. Only one of the 20 final extension sequences contained the shared prefix: the technical discovery branch usually licensed a different extension.

The v7 rule stops at the exact shared 13-base prefix, while retaining the previous one-substitution allowance for the longer cores. It checks both strands, both traversal directions and the target/flank boundary. All 13 bases must match in the shorter motif. The preceding clean prefix can still qualify if it satisfies the existing context requirements; other branches continue. No target IDs, reference sequences or alignment results enter this rule.

The prefix is shared by the published [Illumina adapter sequences](https://support.illumina.com/ko-kr/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html). Recognizing it does not require knowing a sample's downstream adapter variant.

All 21 recorded context routes were traced through the actual graph, verifying 15,343 directed edges. The repaired barrier stops 20 routes after only 12–22 clean added bases, too little to establish context; the remaining route retains 308 clean bases before its adapter and remains eligible as a shorter context. For example, `Q000012` has count 62 at added position 22, then chooses a count-1 `T` edge at position 23 while the other `G` edge has count 86. The old check declared this path external at position 31. The repaired check recognizes the adapter and retains only the preceding 22 bases. A low count alone is not the rejection criterion.

## Complete ordinary-discovery comparison

Same 530 input spacer observations, 398 distinct canonical queries, k=23 graph, two threads, 750-base maximum flank and unchanged search/scoring settings. Reverse complements are consolidated in the same way in both runs.

| Saved output | v6 total | v6 matches | Repaired total | Repaired matches |
|---|---:|---:|---:|---:|
| Protospacer target regions | 172 | 52 | 159 | 53 |
| Extended regions | 154 | 49 | 142 | 50 |

Targets require at least 90% identity and 85% saved-query coverage; extensions require at least 90% identity and 50% whole saved-query coverage in one HSP. BLAST runs only after prediction against the same local curated database. Extension hits need not cover the target interval, matching the established evaluation convention.

The number of nonmatching target sequences falls from 120 to 106; nonmatching extension sequences fall from 105 to 92. Extension reference agreement rises from 31.8% to 35.2%. This is a development-set comparison, not a validated estimate of biological precision.

At the input-query level, all 61 queries that previously had a matching extension retain a matching extension, and one more query gains one: **one gain, zero losses**. All 56 previously matching target queries remain supported; one additional query gains a matching target.

Thirteen of the 21 adapter-affected queries lose their target call. One switches from its adapter-supported exact placement to a different short placement that matches the reference, but has insufficient extension support. Seven retain other qualifying contexts. This is a graph-search rerun, not deletion of a predetermined list of outputs.

The additional matching extension belongs to query `Q000386`: 703/780 bases align at 100% identity (90.1% coverage). Its target and selected side are unchanged. Removing technical branches changes the bounded candidate pool and permits a lower-cost matching continuation to be selected; the precise old pruning event was not instrumented for this query.

Sequence-level changes and query-level changes have different denominators: all 52 old matching target strings remain; one is added. Of the 49 old matching extension strings, 47 remain identical, two are replaced, and three new strings match. None of their previously supported input queries loses its matching-extension outcome.

All 50 matching extensions are 779–792 bases long. The complete saved-extension set contains 138 sequences of 750–799 bases and four of 314–499 bases. These are extended regions, not counts of complete phage genomes.

## Original supplied targets

All 81 original FASTA records are preserved, representing 75 canonical sequences. The supplied-target FASTA, detected-target FASTA and extension FASTA are byte-for-byte identical to v6. All 75 canonical targets extend, producing 69 unique extensions; 43 extensions match, serving 48 original input records. Matching-input gains and losses are both zero.

These supplied targets were selected historically using reference evidence; they are a preservation cohort, separate from the unselected 398-query discovery test.

## Shuffled controls and integrity

The same 398 dinucleotide-preserving shuffled queries were rerun. Their target and extension FASTAs are byte-for-byte unchanged: eight target regions and eight extensions, with one target match and one extension match belonging to different controls. The 8/398 call rate is a background-control observation, not an estimate of precision or false-discovery rate in the real inputs.

The final v7 binary passes the expanded synthetic fixture and its output-schema checks. Native self-tests, batch-helper tests, Python compilation and the diff check pass. All 27 original graph files and the array report retain their original SHA-256 hashes. Only one real-graph process ran at a time, with two threads and a 3.5-GiB process-memory guard.

## Why a stronger abundance filter was not applied

Among the previous 154 extension outputs, the 105 nonmatches were typically smoother than the 49 matches: median adjacent log-count change energy was 0.00492 versus 0.0280. These summaries use the lowest-query-ID profile as the representative of each unique output. Requiring every scored edge to have multiplicity above one would retain only 14 of the 49 matching extensions. None of the 105 nonmatches had an all-count-one profile. Greater smoothness or a blanket minimum-count rule is therefore not a justified fix for these failures.

Additional sequence audits found 23 accepted queries whose strict context later touches known source repeat/junction k-mers; five have matching targets and thirteen have matching extensions. Rejecting every such contact would also discard supported results. Likewise, partial-prefix source-repeat rejection would remove many matching outcomes. These audits were not turned into production filters.

## Scope and implementation

The verified repair addresses missed technical sequence in discovery and extension. **106 target regions and 92 extensions still do not meet the reference criteria.** They are not all explained or validated by this repair; absence from this database alone does not establish their biological origin.

The implementation remains one native class, with no new predictor dependency or reference access. `adapter_position` is the changed prediction function. `write_outputs` now includes `query_id` in `coverage_evidence.tsv`, distinguishing multiple queries that map to one deduplicated target. The native and batch implementation markers are `protospacer-workflow-v7`.

The full-cohort runs froze the repaired prediction code before the final version-marker and coverage-table-label changes. A source comparison verifies that those output-only changes are the only differences from final v7 prediction code. The final binary is also tested on the synthetic fixture.

The regression fixture includes the original one-substitution adapter and a divergent downstream adapter with the exact common prefix. On the identical fixture graph, the previous production binary reports one false new discovery and the repaired binary reports zero. Disabling artifact filtering exposes the expected false context. The fixture also checks preserved supplied records, short-target coordinates, one source passage, source re-entry, pure-spacer membership and a selected single flank. Two intentionally ambiguous repeated-source examples still permit graph-valid routes that are not substrings of the original synthetic molecules; the fixture records this limitation rather than counting them as correct reconstructions.

Local discovery time was 536 seconds versus 472 seconds for v6, including graph loading and index verification and excluding BLAST; search/filter time was 434 versus 390 seconds. These single runs are not a controlled performance benchmark. The existing graph is read-only, and no raw reads were downloaded or rebuilt for this comparison.
