# Reference windows for the SRR14017553 diagnostic audit

This fixture contains all 1,017 qualifying reference HSP locations for the frozen
81 input sequences, plus 24 additional locations where a previously matching
extension HSP spans the complete seed. These latter sites matter because a long
alignment can qualify even when the short seed alone falls below its threshold.
A fresh BLAST of all 81 inputs reproduced the 1,017 seed HSPs exactly.

Seed searches used the local Shah/EarlyVir database (10,021 sequences), BLASTN
`blastn-short`, both strands, word size 7, no DUST/soft masking, E-value ceiling
1000, >=90% identity and >=85% single-HSP query coverage. The additional extension
HSPs use >=90% identity and >=50% single-HSP coverage. No HSP pooling is used.
Checksums and the graph k are recorded in `manifest.json`.

`reference_cases.tsv.gz` contains oriented reference windows up to 750 bases on
each side of each seed location. Coordinates are zero-based, half-open in the
oriented reference sequence. Unaligned query ends are projected from the HSP
boundaries; this projection is a diagnostic hypothesis, not a new alignment.
`reference_allele` is the actual reference interval. `left_outward` is reversed
so its first base is adjacent to the supplied seed; `right_outward` is forward.
Ambiguous reference bases terminate the known flank rather than becoming graph
gaps. Reference-end clipping is explicitly marked. The graph search starts
from the immutable original input endpoint, never a substituted reference allele.

These sequences are used only by separate offline audit executables. MCAAT's
production search and ranker do not load this fixture or the reference database.
The oracle's witnesses prove connected graph paths, not biological provenance.
Its positive paths can include visits or lengths outside production limits;
production-policy exclusions must be reported separately.

Example (from the repository root, after the normal MCAAT build is configured):

```sh
mkdir -p audit
cp bench/fixtures/protospacer-reference-audit/targets.tsv audit/
gzip -dc bench/fixtures/protospacer-reference-audit/reference_cases.tsv.gz > audit/reference_cases.tsv
cmake -S . -B build -DBUILD_PROTO_AUDIT=ON
cmake --build build --target mcaat_reference_audit mcaat_score_reference_flanks --parallel 4
./build/mcaat_reference_audit /path/to/existing/graph/graph audit 1000000
```

The reference oracle emits exact-edge failures, verified path sequences/node IDs,
and bounded edit-search outcomes. `no_path_at_90pct_anchored_identity` is emitted
only with a valid cost lower bound; capped searches without such a certificate
are `search_limit_unresolved`. A million states is a work cap, not proof that
no larger approximate path exists. Exact reference traversal is deterministic
and is checked independently of this cap.

For production decision traces, run
`python3 bench/prepare_protospacer_loss_audit.py audit`, then configure with
`-DPROTO_AUDIT_INSTRUMENTED_INCLUDE=/absolute/path/to/audit/instrumented` and build
`mcaat_trace_extensions`. The generated headers add logging without modifying
the production files. Always verify its selected sequences/coordinates against
an ordinary MCAAT run before attributing losses from the trace.
