#!/usr/bin/env bash
# Build the branch-decision training table from the reference-confirmed
# protospacer regions, then say whether anything in it separates the classes.
#
# For every sample the label file names, label_features maps that sample's
# labelled regions into its graph, walks each region outward from the spacer,
# and writes one row per neighbour at every node with more than one neighbour:
# target=1 for the neighbour the reference continues into, target=0 for its
# siblings. The graph path is read per sample from evidence19_seeds/
# run_summary.json, as it differs between samples.
#
# The label file is not in the repository. Upload it and point LABELS at it.
#
# Writes into OUT:
#   features/<sample>.tsv               one table per sample
#   features/<sample>.tsv.summary.json  that run's counts
#   features/<sample>.log
#   features_all.tsv                    every sample, one header
#   analysis.txt                        class balance, per sample, per feature
#
# Usage:
#   LF=build/label_features LABELS=/vol/d/data/copsac/labels_by_dataset.json \
#       bash scripts/run_label_features.sh                  # every sample
#   LF=... LABELS=... bash scripts/run_label_features.sh SRR14017565
#
# Settings (environment):
#   ROOT=/vol/d/data/copsac  SEEDS=evidence19_seeds  OUT=$ROOT/label_features
#   EXTRA="--lookahead 500"  FORCE=1
set -uo pipefail

LF=${LF:?set LF to the built label_features}
LABELS=${LABELS:?set LABELS to labels_by_dataset.json}
ROOT=${ROOT:-/vol/d/data/copsac}
SEEDS=${SEEDS:-evidence19_seeds}
OUT=${OUT:-$ROOT/label_features}
EXTRA=${EXTRA:-}
FORCE=${FORCE:-0}

[[ -x $LF || -f $LF ]] || { echo "no label_features at $LF" >&2; exit 2; }
[[ -f $LABELS ]] || { echo "no label file at $LABELS" >&2; exit 2; }

field() { python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get(sys.argv[2], ""))' "$1" "$2" 2>/dev/null; }

mkdir -p "$OUT/features"

# ── 1. which samples ─────────────────────────────────────────────────────────
if (( $# )); then
  samples=("$@")
else
  mapfile -t samples < <("$LF" --labels "$LABELS" --list | awk -F'\t' 'NF==2 && $2+0>0 {print $1}')
fi
(( ${#samples[@]} )) || { echo "no samples to run" >&2; exit 1; }
echo "label file $LABELS"
"$LF" --labels "$LABELS" --list | head -1
echo "samples to run: ${#samples[@]}"
echo

# ── 2. one run per sample ────────────────────────────────────────────────────
ran=0 skipped=0
for s in "${samples[@]}"; do
  s=${s%/}; s=$(basename "$s")
  table=$OUT/features/$s.tsv
  if [[ -s $table && $FORCE != 1 ]]; then
    echo "$s: already done, skipped"; (( skipped++ )); continue
  fi
  # The graph is not in one place across samples: most have it under
  # protospacer_batch/graph, a few have it at graph/, and some have none at
  # all. The candidates are tried in order and the one used is always named,
  # because mapping a sample's labels into another sample's graph would not
  # announce itself in the output.
  summary=$ROOT/$s/$SEEDS/run_summary.json
  graph=""; source=""
  declared=$(field "$summary" graph)
  if [[ -n $declared && -f ${declared}.sdbg_info ]]; then
    graph=$declared; source=$SEEDS/run_summary.json
  else
    [[ -n $declared ]] && echo "$s: $SEEDS/run_summary.json names $declared, which is not there"
    for candidate in "$ROOT/$s/protospacer_batch/graph/graph" "$ROOT/$s/graph/graph"; do
      [[ -f ${candidate}.sdbg_info ]] || continue
      graph=$candidate; source="found at ${candidate#"$ROOT/$s/"}"
      break
    done
  fi
  if [[ -z $graph ]]; then
    echo "$s: no graph found (tried $SEEDS/run_summary.json, protospacer_batch/graph, graph), skipped"
    (( skipped++ )); continue
  fi
  echo "=== $s  graph=$graph  [$source]"
  # shellcheck disable=SC2086
  "$LF" --graph "$graph" --labels "$LABELS" --dataset "$s" --output "$table" $EXTRA \
      2>&1 | tee "$OUT/features/$s.log"
  if [[ -s $table ]]; then (( ran++ )); else echo "$s: no rows written"; fi
done
echo
echo "ran $ran, skipped $skipped"

# ── 3. one table ─────────────────────────────────────────────────────────────
all=$OUT/features_all.tsv
first=1
: > "$all"
for s in "${samples[@]}"; do
  s=$(basename "${s%/}")
  table=$OUT/features/$s.tsv
  [[ -s $table ]] || continue
  if (( first )); then cat "$table" >> "$all"; first=0
  else tail -n +2 "$table" >> "$all"; fi
done
rows=$(( $(wc -l < "$all") - 1 ))
(( rows > 0 )) || { echo "no rows in $all, nothing to analyse" >&2; exit 1; }
echo "combined $rows rows into $all"

# ── 4. analysis ──────────────────────────────────────────────────────────────
# Column 9 is the target. Columns 10 to 18 are the nine features, every one of
# them read from topology and multiplicity. For each, the class means and the
# Mann-Whitney AUC: 0.5 is no separation, and a value far from it on a single
# feature says a rule would do what a model is being considered for.
report=$OUT/analysis.txt
{
  echo "label_features analysis"
  echo "table   $all"
  echo "labels  $LABELS"
  echo "written $(date -u '+%Y-%m-%d %H:%M:%SZ')"
  echo

  echo "── run counts, from the per-sample summaries ──"
  python3 - "$OUT/features" <<'PY'
import json, os, sys
keys = ["entries", "entries_mapped", "sides_mapped", "sides_empty", "anchors",
        "reanchors", "steps", "decisions", "decisions_skipped", "positives",
        "negatives", "rows"]
total = dict.fromkeys(keys, 0)
folder = sys.argv[1]
for name in sorted(os.listdir(folder)):
    if not name.endswith(".summary.json"):
        continue
    with open(os.path.join(folder, name)) as handle:
        d = json.load(handle)
    for k in keys:
        total[k] += int(d.get(k, 0) or 0)
width = max(len(k) for k in keys)
for k in keys:
    print("  %-*s %12d" % (width, k, total[k]))
if total["decisions"]:
    print("  %-*s %12.1f%%" % (width, "decisions skipped",
                               100.0 * total["decisions_skipped"] / total["decisions"]))
if total["entries"]:
    print("  %-*s %12.1f%%" % (width, "regions mapped",
                               100.0 * total["entries_mapped"] / total["entries"]))
PY
  echo

  echo "── class balance ──"
  awk -F'\t' 'NR>1{n++; if($9==1) p++; else q++}
       END{printf "  rows %d, positives %d (%.2f%%), negatives %d, negatives per positive %.2f\n",
                  n, p, 100*p/n, q, (p?q/p:0)}' "$all"
  echo

  echo "── decisions and rows per sample ──"
  printf "  %-14s %9s %10s %10s %9s %9s\n" sample rows positives negatives labels spacers
  awk -F'\t' 'NR>1{r[$1]++; if($9==1) p[$1]++; else q[$1]++;
                   lab[$1"\t"$2]=1; sp[$1"\t"$3]=1}
       END{for (s in r) {
             nl=0; ns=0;
             for (k in lab) { split(k,a,"\t"); if (a[1]==s) nl++ }
             for (k in sp)  { split(k,a,"\t"); if (a[1]==s) ns++ }
             printf "  %-14s %9d %10d %10d %9d %9d\n", s, r[s], p[s]+0, q[s]+0, nl, ns }}' \
      "$all" | sort
  echo

  echo "── per feature: class means and single-feature AUC ──"
  echo "   AUC is the chance a positive neighbour outranks a negative one at the"
  echo "   same decision set. 0.50 is nothing; far from it means a rule suffices."
  echo
  header=$(head -1 "$all")
  printf "  %-22s %14s %14s %7s %8s\n" feature mean_positive mean_negative AUC "|AUC-.5|"
  for col in $(seq 10 18); do
    name=$(echo "$header" | cut -f"$col")
    means=$(awk -F'\t' -v c="$col" 'NR>1{ if($9==1){sp+=$c; np++} else {sq+=$c; nq++} }
            END{printf "%.4f %.4f", (np?sp/np:0), (nq?sq/nq:0)}' "$all")
    auc=$(tail -n +2 "$all" | cut -f9,"$col" | sort -k2,2g | awk '
            {t[NR]=$1; v[NR]=$2}
            END{n=NR; i=1;
                while (i<=n) { j=i; while (j<n && v[j+1]==v[i]) j++;
                               r=(i+j)/2;
                               for (m=i;m<=j;m++) { if (t[m]==1){s+=r; np++} else nq++ }
                               i=j+1 }
                if (np==0||nq==0) { printf "nan" } else { printf "%.4f", (s-np*(np+1)/2)/(np*nq) }}')
    # shellcheck disable=SC2086
    printf "  %-22s %14s %14s %7s %8s\n" "$name" $means "$auc" \
        "$(awk -v a="$auc" 'BEGIN{ if (a=="nan") print "nan"; else printf "%.4f", (a>0.5?a-0.5:0.5-a) }')"
  done | sort -k5,5gr
  echo
  echo "Splitting train from test: split on spacer_key, or on dataset, never on"
  echo "rows. One spacer contributes every neighbour at every decision on its"
  echo "region, so a row-wise split puts the same decision on both sides."
} > "$report" 2>&1

echo "wrote $report"
echo
sed -n '/── class balance/,$p' "$report" | head -60

# Fitting a model is not done here: the table is the deliverable, and the
# analysis above says whether a model is worth fitting at all.
