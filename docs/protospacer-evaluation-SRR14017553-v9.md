# SRR14017553: final protospacers require extensions (v9)

V9 changes final output eligibility: a protospacer is saved as a final result only when its selected extension passes the existing filters. It does not change discovery, path selection, extension thresholds or the coverage filter introduced in v8.

The full native workflow was rerun on the same existing k=23 graph: 49 arrays, 530 source-spacer observations and 398 distinct corrected spacers, with two threads. No reference search was run for this output-policy check.

| Output | V8 | V9 |
|---|---:|---:|
| All candidate protospacer regions | 159 | 159 |
| Final `Protospacers.fasta` | 159 | 117 |
| Candidate regions without an accepted extension | 42 | 42 |
| Accepted `ContiguousPaths.fasta` extensions | 106 | 106 |
| Queries with an accepted extension | 119 | 119 |

All 106 accepted extension sequences and IDs are byte-identical to v8. The 117 final target sequences are exactly the subset linked to accepted extensions. Multiple queries can share a target, and multiple targets can share an extension; these are different counting units.

The 42 candidate-only regions remain in the diagnostic output. This restriction narrows the final results by extension support; it does not establish that every retained region is biologically correct.

## Output contract

- `Protospacers.fasta`: final target regions, each linked to at least one saved extension.
- `ContiguousPaths.fasta`: accepted extensions; the existing 200–750-base flank limits and other filters still apply.
- `associations.tsv`: final source/query/target/extension links with populated extension IDs and exact coordinates.
- `CandidateProtospacers.fasta` and `candidate_associations.tsv`: every detected candidate and its associations, including final results. These are diagnostic files; their counts must not be added to final counts.
- `SuppliedProtospacers.fasta`: every supplied input record, preserved unchanged even if it has no accepted extension.

Target IDs are shared between diagnostic and final files; gaps in final numbering are expected.

## Verification

The final-target subset, every association coordinate and all summary counts pass validation. Sources, query outcomes, search work, selected graph routes, coverage evidence, withheld routes and original candidate associations remain unchanged. All 28 original graph/input component hashes are unchanged.

Eight native fixture outputs pass validation, including source-array passage/re-entry, short targets and reverse complements, exact flank limits, and the all-extensions-withheld case. Supplied input preservation is checked by the source-context fixture. A batch fixture verifies the completed `NO_ACCEPTED_EXTENSIONS` outcome and resumption without rerunning it. The batch validator also checks that multiple targets may share one extension and rejects a final target missing its extension association.

No additional biological accuracy analysis is claimed for this output-policy change. See the [v8 reference evaluation](protospacer-evaluation-SRR14017553-v8.md) for the unchanged extension sequences' previous results, and the [workflow guide](protospacer-extension.md) for the server command and function descriptions.
