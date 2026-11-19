#!/usr/bin/env bash
# Re-run discovery and the 0.998 support stage for every sample whose original
# evidence19 stages completed, using the same parameters as those runs. Only the
# new binary's defaults differ (source-array crossing off, CRISPR repeat barrier on).
#
# Writes only new folders per sample (existing folders are read, never changed):
#   SRR*/<LABEL>_seeds/  SRR*/<LABEL>_support_0.998/  and their .log files
#
# Usage:
#   MCAAT=/path/to/build/mcaat bash run_evidence19_stages.sh                # all completed samples
#   MCAAT=/path/to/build/mcaat bash run_evidence19_stages.sh SRR14017553    # named samples only
# Settings (environment): ROOT=/vol/d/data/copsac THREADS=28 LABEL=v10
#   BARRIER=true|false  (false = crossing off only, for the ablation run; use LABEL=v10nobarrier)
# Resumable: completed stages are skipped; an incomplete new folder is renamed
# to *.incomplete.<time> and that stage is run again.
set -uo pipefail
MCAAT=${MCAAT:?set MCAAT to the new build/mcaat}
ROOT=${ROOT:-/vol/d/data/copsac}
THREADS=${THREADS:-28}
LABEL=${LABEL:-v10}
BARRIER=${BARRIER:-true}

field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2" 2>/dev/null; }
complete() { [[ -f $1 && $(field "$1" status) == complete ]]; }
stamp() { date '+%F %T'; }

"$MCAAT" phage-curate --help 2>&1 | grep -q -- '--crispr-repeat-barrier' \
  || { echo "$MCAAT is not the target-region build (no --crispr-repeat-barrier); stop."; exit 1; }
"$MCAAT" phage-curate --self-test || { echo "self-test failed; stop."; exit 1; }

if (( $# )); then dirs=(); for s in "$@"; do dirs+=("$ROOT/$s"); done
else dirs=("$ROOT"/SRR*); fi

for dir in "${dirs[@]}"; do
  dir=${dir%/}; s=$(basename "$dir")
  complete "$dir/evidence19_seeds/run_summary.json" || { echo "$s: no complete evidence19_seeds, skipped"; continue; }
  complete "$dir/evidence19_support_0.998/summary.json" || { echo "$s: no complete evidence19_support_0.998, skipped"; continue; }
  graph=$(field "$dir/evidence19_seeds/run_summary.json" graph)
  [[ -n $graph ]] || { echo "$s: graph path missing in run_summary.json, skipped"; continue; }
  seeds=$dir/${LABEL}_seeds
  support=$dir/${LABEL}_support_0.998

  if ! complete "$seeds/run_summary.json"; then
    [[ -e $seeds ]] && mv "$seeds" "$seeds.incomplete.$(date +%s)"
    echo "$(stamp) $s discovery"
    "$MCAAT" phage-curate --arrays "$dir" --graph "$graph" --output "$seeds" --threads "$THREADS" \
      --source-array-crossing false --crispr-repeat-barrier "$BARRIER" \
      > "$dir/${LABEL}_seeds.log" 2>&1 \
      || { echo "$(stamp) $s discovery FAILED, see $dir/${LABEL}_seeds.log"; continue; }
  fi

  if ! complete "$support/summary.json"; then
    [[ -e $support ]] && mv "$support" "$support.incomplete.$(date +%s)"
    echo "$(stamp) $s support 0.998"
    "$MCAAT" phage-curate --arrays "$dir" --graph "$graph" --extension-seeds "$seeds" \
      --output "$support" --coverage-method evidence --min-contig-length 500 \
      --max-flank-length 750 --beam-width 64 --max-context-paths 256 \
      --min-branch-support 0.998 --max-reported-contigs 0 --threads "$THREADS" \
      --source-array-crossing false --crispr-repeat-barrier "$BARRIER" \
      > "$dir/${LABEL}_support_0.998.log" 2>&1 \
      || { echo "$(stamp) $s support FAILED, see $dir/${LABEL}_support_0.998.log"; continue; }
  fi
  echo "$(stamp) $s done"
done
echo "$(stamp) batch finished"
