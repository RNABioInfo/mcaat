#!/usr/bin/env bash
set -euo pipefail

# Simple benchmark runner for MCAAT cycle_finder
# Produces CSV rows: threads,run,elapsed_seconds,max_rss_kb,cycles

usage() {
  cat <<EOF
Usage: $0 --binary PATH --input PATH --threads LIST --runs N -o OUTCSV

Options:
  --binary PATH     Path to cycle_finder binary
  --input PATH      Input graph file
  --threads LIST    Comma-separated thread counts (e.g. 1,4,8,24)
  --runs N          Number of runs per thread count (default: 3)
  -o OUTCSV         Output CSV file (overwritten)

Example:
  $0 --binary ./bin/cycle_finder --input graphs/huge_graph.bin --threads 1,4,8,24 --runs 3 -o bench/results.csv
EOF
  exit 1
}

BINARY="" INPUT="" THREADS="" RUNS=3 OUTCSV=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --binary) BINARY="$2"; shift 2;;
    --input) INPUT="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --runs) RUNS="$2"; shift 2;;
    -o) OUTCSV="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

if [[ -z "$BINARY" || -z "$INPUT" || -z "$THREADS" || -z "$OUTCSV" ]]; then
  usage
fi

mkdir -p "$(dirname "$OUTCSV")"

echo "threads,run,elapsed_s,max_rss_kb,cycles" > "$OUTCSV"

IFS=',' read -r -a THREAD_ARR <<< "$THREADS"

for threads in "${THREAD_ARR[@]}"; do
  for ((run=1; run<=RUNS; run++)); do
    outtmp="/tmp/cycle_out.${threads}.${run}.json"
    timetmp="/tmp/cycle_time.${threads}.${run}.txt"

    echo "Running threads=$threads run=$run"

    # Use /usr/bin/time to capture wall time and max RSS
    /usr/bin/time -f "%e %M" -o "$timetmp" "$BINARY" --input "$INPUT" --threads "$threads" --out "$outtmp"

    read -r elapsed rss_kb < "$timetmp"

    # Try to extract cycles from output JSON if possible
    cycles=""
    if command -v jq >/dev/null 2>&1 && jq -e . "$outtmp" >/dev/null 2>&1; then
      # Common fields: either cycles array or cycle_count
      if jq -e '.cycles' "$outtmp" >/dev/null 2>&1; then
        cycles=$(jq '.cycles | length' "$outtmp")
      elif jq -e '.cycle_count' "$outtmp" >/dev/null 2>&1; then
        cycles=$(jq '.cycle_count' "$outtmp")
      fi
    else
      # Fallback: grep for numbers labeled cycles or count keys
      cycles=$(grep -Eo '"cycle_count"[[:space:]]*:[[:space:]]*[0-9]+' "$outtmp" | head -n1 | grep -Eo '[0-9]+') || true
      if [[ -z "$cycles" ]]; then
        cycles=$(grep -Eo '"cycles"[[:space:]]*:[[:space:]]*[0-9]+' "$outtmp" | head -n1 | grep -Eo '[0-9]+') || true
      fi
    fi

    echo "${threads},${run},${elapsed},${rss_kb},${cycles}" >> "$OUTCSV"

    rm -f "$outtmp" "$timetmp"
  done
done

echo "Done. Results at: $OUTCSV"