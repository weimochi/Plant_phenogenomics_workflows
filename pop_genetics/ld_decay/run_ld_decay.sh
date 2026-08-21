#!/usr/bin/env bash
set -euo pipefail
[[ $# -eq 5 ]] || { echo "Usage: $0 BFILE GROUP_DIR_OR_DASH OUT_DIR MAX_DISTANCE_KB THREADS" >&2; exit 2; }
BFILE=$1; GROUP_DIR=$2; OUT_DIR=$3; MAX_KB=$4; THREADS=$5
command -v plink >/dev/null || { echo "ERROR: plink not found" >&2; exit 2; }
for ext in bed bim fam; do [[ -s "${BFILE}.${ext}" ]] || { echo "ERROR: missing ${BFILE}.${ext}" >&2; exit 2; }; done
[[ "$MAX_KB" =~ ^[1-9][0-9]*$ && "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: distance and threads must be positive integers" >&2; exit 2; }
mkdir -p "$OUT_DIR"

run_ld() {
    local name=$1; shift
    echo "[$(date)] LD: $name"
    plink --bfile "$BFILE" --allow-extra-chr "$@" --r2 gz \
        --ld-window 999999 --ld-window-kb "$MAX_KB" --ld-window-r2 0 \
        --threads "$THREADS" --out "$OUT_DIR/$name"
    [[ -s "$OUT_DIR/$name.ld.gz" ]] || { echo "ERROR: missing LD output for $name" >&2; exit 3; }
}

if [[ "$GROUP_DIR" == "-" ]]; then run_ld whole_cohort; exit 0; fi
[[ -d "$GROUP_DIR" ]] || { echo "ERROR: group directory not found" >&2; exit 2; }
mapfile -t GROUP_FILES < <(find "$GROUP_DIR" -maxdepth 1 -type f -name '*.keep' | sort)
(( ${#GROUP_FILES[@]} >= 1 )) || { echo "ERROR: no *.keep files found" >&2; exit 2; }
for group in "${GROUP_FILES[@]}"; do
    [[ -s "$group" ]] || { echo "ERROR: empty group: $group" >&2; exit 2; }
    [[ $(awk 'NF != 2 {bad++} END {print bad+0}' "$group") -eq 0 ]] || { echo "ERROR: keep file must have FID IID: $group" >&2; exit 2; }
    run_ld "$(basename "$group" .keep)" --keep "$group"
done
