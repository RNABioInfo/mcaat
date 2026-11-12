# SRR14017553: all 81 inputs, target placement, and extension recovery

**Superseded by the [complete-placement audit](../protospacer-complete-placement-2026-11-12/README.md): 12 demonstrated recoverable inputs, 69 completed negative classifications, and 11 production matches. The results below describe the earlier run.**

The earlier claim that almost all 81 extensions could be recovered was wrong.
Checking the starting sequence, the intervening edges, and the emitted DNA
shows **11 inputs with a demonstrated reference-consistent continuation**.
The updated graph-only software selects a matching extension for **10**, up
from **8**. **One selection remains unresolved.** Eleven is a demonstrated
count, not a proved maximum over every possible strain or approximate path.

All 81 input records remain byte-for-byte unchanged. Counts below refer to
input records and their extensions, not distinct phages or independent genomes.

## Measured production result

The existing k=23 graph was reused. Each final beam-width run used 24 threads,
750 bp maximum per flank, 1,550 bp maximum total, and at most three visits to a
node outside source-array cycles. References were used only after prediction
and in separate diagnostic executables.

| Run | Inputs retained | Inputs extended | Unique output regions | Matching extensions | Previous matches lost |
|---|---:|---:|---:|---:|---:|
| Previous N=32 | 81 | 45 | 40 | 8 | — |
| New N=16 | 81 | 54 | 51 | 10 | 0 |
| New N=32 | 81 | 56 | 51 | 10 | 0 |
| New N=64 | 81 | 56 | 53 | 10 | 0 |

The local BLAST criterion is >=90% identity and >=50% single-HSP coverage of
the emitted extension. The two additional matches are:

- `seed_7a5827b9089bd147`: **1,537 bp**, matching OTU_2191 across the entire
  sequence at 100% identity.
- `seed_c680dd8ddedb88b2`: **133 bp**, matching OTU_12556 across the entire
  sequence at 100% identity. This is a short partial extension, not a complete phage.

The other output regions remain predictions; a missing database match does not
by itself establish that they are false positives.

## What failed and what changed

The supplied spacer and the actual target allele can differ, particularly at
unaligned query ends. The earlier retained-input mode required the entire
original string to be an exact graph path. It consequently anchored some queries
in their source array and could not enter a matching viral allele.

The implementation now keeps any existing exact-placement extension unchanged.
If there is none, it searches actual graph sequences for a related target at
>=90% identity over >=85% of the query. A short query can be placed within a
real graph k-mer using BOSS suffix lookup. Every surrounding base comes from
observed graph labels/edges. The original FASTA is preserved; the observed allele
and its aligned query coordinates are reported separately.

All candidate placements first receive a short flank probe. The best N then
receive the remaining extension budget. This prevents the first repetitive
placement from exhausting the budget before other locations are considered.
Source-array, adapter, cycle, multiplicity and length guards still apply.
The mapper is bounded and uses substitutions with exact/one-substitution anchors;
it is not a complete indel-aware target mapper.

The eleven graph-produced matching sequences are in
[reference_matching_graph_extensions.fasta](reference_matching_graph_extensions.fasta).
This is a **diagnostic collection chosen by reference evaluation**: ten are final
N=32 predictions; the eleventh comes from the isolated-placement experiment. It
is not the ordinary production output.

## The one remaining selection failure

`seed_12839d3c8ccafc84` is a 17-base input. A related target at OTU_12626 is present
in the graph. When its independently discovered graph placement is extended in
isolation, N=32 produces **1,070 bp matching that reference exactly**: 303 bp
left, the 17-base observed target, and 750 bp right. The left side then reaches
an adapter-containing continuation and is stopped.

The ordinary run considers many other matching locations. The correct location
survives the placement/probe stage, but another location has a better multiplicity
score and is selected. That selected location does not match the local reference.
Increasing N to 64 does not resolve this choice. The measured median edge counts
are about **230–232** on the selected path, with coefficients of variation
**0.04–0.05**. On the reference-consistent path the medians are **3–5**, with
coefficients of variation **0.35–0.45**. The current uniformity score therefore
favors the competing location. These values come from the actual graph, recorded
in [remaining_ambiguity_multiplicity.json](remaining_ambiguity_multiplicity.json). The result is explicitly marked
ambiguous; the diagnostic placement table retains the competing sequences and
scores. **The software does not yet select this eleventh known matching extension.**
Forcing the reference-supported location to win would make this benchmark
reference-guided, which is outside the requested graph-only method.

## Every input is accounted for

[all81_results.tsv](all81_results.tsv) gives one row for every original input.
The tested starting placements comprise each original input and every target
allele at all 1,041 reference locations (1,017 qualifying seed HSPs plus 24
previous extension-derived locations).

- **10:** a matching continuation exists and production selects a matching extension.
- **1:** a matching continuation exists, but production selects another target location.
- **33:** at least one exact reference allele maps, but no qualifying continuation
  exists at those tested locations, or from the original input, under the audit rules.
- **37:** the reference alleles do not map; original-input searches also find no
  qualifying continuation under the audit rules.

The continuation test allows substitutions/indels at >=90% aligned identity to
any known reference prefix from 100 through 750 bases. It excludes default
adapters and closed walks wholly within source-array cycles. It does not use the
production multiplicity filter, beam, or initial repeat-similarity filter.
Negative certificates require a completed search; a work cap is never absence.
Short queries are checked across possible graph contexts. Unknown reference bases
and contig ends terminate the known comparison interval. A 100-reference-base
alignment can emit fewer graph bases when it contains deletions.

One stricter search requiring the unchanged original sequence of the already
recovered `7a58...` input remains work-limited; the demonstrated variant path
answers the existence question for that input. Some alternative reference
locations of other positive inputs remain work-limited too. The report does
**not** prove that no untested target variant or unknown reference flank could
extend further. [summary.json](summary.json) records this limitation explicitly.

## Verification and reproduction

The release and ASan/UBSan suites pass. Tests compare the diagnostic algorithms
with independent string DP and exhaustive small graph walks, including backward
walking, indels, source cycles, short seeds and work limits. Native BOSS suffix
results were checked against exhaustive padding enumeration for the five short
inputs. The 2- and 24-thread predictions are identical. The frozen input checksum is
`9a7886832f257d1d7e4bc69f0cd64aab9354195dbdad47a456e5c6b1376cff12`.

[validation.json](validation.json) records versions/checksums; compressed
per-location certificates and path witnesses are included beside this report.
The graph files were not rebuilt or changed. For the complete server command,
see [rerun.sh](rerun.sh).
