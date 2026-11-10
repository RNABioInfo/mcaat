# Observed SRR14017553 graph fragments

These k=23 edges and multiplicities were read from the downloaded SRR14017553
SDBG on 2026-11-10. The SHA-256 of `graph.sdbg_info` was
`94ca5ad4422b97de35f7e0d46e165bdfb8efaaabe283166ecd3ea69360e21bce`.

`cases.tsv` contains two unchanged accepted protospacers, their source-array
contexts, and reference-supported sequences that the original graph contains.
`edges.tsv` contains observed edges along the diagnostic paths and their immediate
alternatives. Counts were not simulated or adjusted. The restricted graph can
have artificial dead ends at the boundary of this extraction, so this is a
regression fixture, not a replacement for evaluation on the complete graph.

Before the fixes, the first supplied interval stopped at a 1,315-base extension.
Suppressing artificial terminal prefixes alone reached 1,450 bases; following
gradual local depth changes as well reaches the complete 1,530-base target.
The second complete target is 1,528 bases. The tests compare the complete emitted
strings, not just lengths, and preserve the supplied interval within each string.

The production algorithm does not read these fixtures or reference sequences.
See `docs/protospacer-extension.md` for the full-graph before/after results and
their limits.
