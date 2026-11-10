# Preserved protospacer baseline

These are the 81 boundary-corrected exported intervals with a local reference
match at >=90% identity and >=85% single-alignment query coverage in the
2026-11-10 SRR14017553 audit. IDs, order and sequences are preserved. The rules
used to propose boundaries did not depend on reference matches. Seventy-one
records also meet E <=0.01. These are reference-matching intervals, not confirmed
phages or a precision estimate for the caller.

The fixture is input data, never a hardcoded exception in MCAAT. In particular,
`--protospacers` must retain any supplied FASTA independently of extension success.
Five records are shorter than k=23: 17, 17, 19, 21 and 22 bases. They are retained
with `shorter_than_graph_k` if extension on a k=23 graph is requested. No missing
bases are invented to make an anchor.

The graph used by `validate_retained_protospacers.py` is deliberately unrelated
to these 81 records. That portion checks preservation under extension failure,
not recovery of their real flanks. Separate synthetic targets in the same suite
verify traversal, sequence truth, adapters, ambiguity and thread determinism.
The original pooled graph must be tested on the server.
