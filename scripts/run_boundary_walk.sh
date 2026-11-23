#!/usr/bin/env bash
# Run boundary_walk over the evidence19 seed sources of one or more samples.
#
# The graph path is read per sample from evidence19_seeds/run_summary.json, as
# it differs between samples. Every setting is left at the binary's defaults,
# which are the configuration measured on SRR14017565:
#   --window-low 0.1 --window-high 5 --window-running true
#   --min-harmonic-mean 3 --require-connected-spacer false
# Pass EXTRA to override any of them for an ablation.
#
# Writes per sample:  SRR*/<LABEL>/<criterion>.fasta  summary.tsv  spacers.tsv
# and SRR*/<LABEL>.log. An existing folder is left alone unless FORCE=1.
#
# Usage:
#   BW=/path/to/build/boundary_walk bash run_boundary_walk.sh SRR14017570
#   BW=... bash run_boundary_walk.sh                       # every sample with seeds
# Settings (environment): ROOT=/vol/d/data/copsac LABEL=bw SEEDS=evidence19_seeds
#   EXTRA="--min-harmonic-mean 5"   FORCE=1
set -uo pipefail
BW=${BW:?set BW to the built boundary_walk}
ROOT=${ROOT:-/vol/d/data/copsac}
LABEL=${LABEL:-bw}
SEEDS=${SEEDS:-evidence19_seeds}
EXTRA=${EXTRA:-}
FORCE=${FORCE:-0}

field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2" 2>/dev/null; }

if (( $# )); then dirs=(); for s in "$@"; do dirs+=("$ROOT/$s"); done
else dirs=("$ROOT"/SRR*); fi

for dir in "${dirs[@]}"; do
  dir=${dir%/}; s=$(basename "$dir")
  sources=$dir/$SEEDS/sources.tsv
  summary=$dir/$SEEDS/run_summary.json
  out=$dir/$LABEL
  [[ -f $sources ]] || { echo "$s: no $SEEDS/sources.tsv, skipped"; continue; }
  graph=$(field "$summary" graph)
  [[ -n $graph ]] || { echo "$s: no graph path in $SEEDS/run_summary.json, skipped"; continue; }
  if [[ -s $out/summary.tsv && $FORCE != 1 ]]; then echo "$s: $LABEL already done, skipped"; continue; fi
  (( FORCE == 1 )) && rm -rf "$out"
  echo "=== $s  graph=$graph"
  # shellcheck disable=SC2086
  "$BW" "$graph" "$sources" "$out" $EXTRA 2>&1 | tee "$dir/$LABEL.log"
  awk -F'\t' 'NR>1 && $1=="harmonic_mean"{n++} END{printf "   %s: %d contigs kept\n", "'"$s"'", n+0}' \
      "$out/summary.tsv" 2>/dev/null
done

cat <<'NEXT'

Scoring runs on the machine with the reference database, not here:
  python3 eval_criteria.py --dir SRR*/LABEL --db all_curated_db --out EVAL
  python3 rule_table.py --dir SRR*/LABEL --eval EVAL --criterion harmonic_mean
NEXT
