#!/usr/bin/env bash
set -euo pipefail
cd /vol/d/mcaat-protospacer/mcaat-protospacer-full
git switch codex/protospacer-extension
git pull --ff-only origin codex/protospacer-extension
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 24
./build/mcaat phage-curate \
  --arrays /vol/d/data/copsac/SRR14017553 \
  --graph /vol/d/data/copsac/SRR14017553/graph/graph \
  --protospacers bench/fixtures/protospacer-baseline-81/Protospacers.fasta \
  --output /vol/d/data/copsac/SRR14017553/protospacer_variant_extensions \
  --threads 24 --beam-width 32 --max-node-visits 3 \
  --max-flank-length 750 --max-length 1550 --max-states 50000
