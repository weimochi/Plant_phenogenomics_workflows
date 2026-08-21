#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 5 ]] || { echo "Usage: $0 VCF GROUP_DIR_OR_DASH OUT_DIR WINDOW_SIZE WINDOW_STEP" >&2; exit 2; }
VCF=$1; GROUP_DIR=$2; OUT_DIR=$3; WINDOW_SIZE=$4; WINDOW_STEP=$5
command -v vcftools >/dev/null || { echo "ERROR: vcftools not found" >&2; exit 2; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found" >&2; exit 2; }
[[ -s "$VCF" ]] || { echo "ERROR: VCF not found: $VCF" >&2; exit 2; }
[[ "$WINDOW_SIZE" =~ ^[1-9][0-9]*$ && "$WINDOW_STEP" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: window parameters must be positive integers" >&2; exit 2; }
mkdir -p "$OUT_DIR"

if [[ "$GROUP_DIR" == "-" ]]; then
    vcftools --gzvcf "$VCF" --window-pi "$WINDOW_SIZE" --window-pi-step "$WINDOW_STEP" --out "$OUT_DIR/whole_cohort"
    [[ -s "$OUT_DIR/whole_cohort.windowed.pi" ]] || { echo "ERROR: missing whole-cohort output" >&2; exit 3; }
    exit 0
fi

[[ -d "$GROUP_DIR" ]] || { echo "ERROR: group directory not found: $GROUP_DIR" >&2; exit 2; }
mapfile -t GROUP_FILES < <(find "$GROUP_DIR" -maxdepth 1 -type f -name '*.samples.txt' ! -name 'admixed.samples.txt' | sort)
(( ${#GROUP_FILES[@]} >= 1 )) || { echo "ERROR: no *.samples.txt files found" >&2; exit 2; }
VCF_SAMPLES=$(mktemp)
trap 'rm -f "$VCF_SAMPLES"' EXIT
bcftools query -l "$VCF" | sort -u > "$VCF_SAMPLES"

for group in "${GROUP_FILES[@]}"; do
    [[ -s "$group" ]] || { echo "ERROR: empty group: $group" >&2; exit 2; }
    missing=$(comm -23 <(sort -u "$group") "$VCF_SAMPLES")
    [[ -z "$missing" ]] || { echo "ERROR: samples absent from VCF in $group: $missing" >&2; exit 2; }
    name=$(basename "$group" .samples.txt)
    echo "[$(date)] Nucleotide diversity: $name"
    vcftools --gzvcf "$VCF" --keep "$group" --window-pi "$WINDOW_SIZE" --window-pi-step "$WINDOW_STEP" --out "$OUT_DIR/$name"
    [[ -s "$OUT_DIR/$name.windowed.pi" ]] || { echo "ERROR: missing output for $name" >&2; exit 3; }
done
