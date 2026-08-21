#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 5 ]] || { echo "Usage: $0 VCF GROUP_DIR OUT_DIR WINDOW_SIZE WINDOW_STEP" >&2; exit 2; }
VCF=$1
GROUP_DIR=$2
OUT_DIR=$3
WINDOW_SIZE=$4
WINDOW_STEP=$5

command -v vcftools >/dev/null || { echo "ERROR: vcftools not found" >&2; exit 2; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found" >&2; exit 2; }
[[ -s "$VCF" ]] || { echo "ERROR: VCF not found: $VCF" >&2; exit 2; }
[[ -d "$GROUP_DIR" ]] || { echo "ERROR: group directory not found: $GROUP_DIR" >&2; exit 2; }
[[ "$WINDOW_SIZE" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid window size" >&2; exit 2; }
[[ "$WINDOW_STEP" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid window step" >&2; exit 2; }
mkdir -p "$OUT_DIR"

mapfile -t GROUP_FILES < <(find "$GROUP_DIR" -maxdepth 1 -type f -name '*.samples.txt' ! -name 'admixed.samples.txt' | sort)
(( ${#GROUP_FILES[@]} >= 2 )) || { echo "ERROR: at least two population lists are required" >&2; exit 2; }

VCF_SAMPLES=$(mktemp)
trap 'rm -f "$VCF_SAMPLES"' EXIT
bcftools query -l "$VCF" | sort -u > "$VCF_SAMPLES"

for group in "${GROUP_FILES[@]}"; do
    [[ -s "$group" ]] || { echo "ERROR: empty group: $group" >&2; exit 2; }
    [[ $(sort -u "$group" | wc -l) -eq $(wc -l < "$group") ]] || { echo "ERROR: duplicate IDs in $group" >&2; exit 2; }
    missing=$(comm -23 <(sort -u "$group") "$VCF_SAMPLES")
    [[ -z "$missing" ]] || { echo "ERROR: samples absent from VCF in $group: $missing" >&2; exit 2; }
done

for ((i=0; i<${#GROUP_FILES[@]}; i++)); do
    for ((j=i+1; j<${#GROUP_FILES[@]}; j++)); do
        g1=${GROUP_FILES[$i]}
        g2=${GROUP_FILES[$j]}
        overlap=$(comm -12 <(sort -u "$g1") <(sort -u "$g2"))
        [[ -z "$overlap" ]] || { echo "ERROR: overlapping groups: $g1 and $g2" >&2; exit 2; }
        n1=$(basename "$g1" .samples.txt)
        n2=$(basename "$g2" .samples.txt)
        out="$OUT_DIR/${n1}_vs_${n2}"
        echo "[$(date)] FST: $n1 vs $n2"
        vcftools --gzvcf "$VCF" --weir-fst-pop "$g1" --weir-fst-pop "$g2" \
            --fst-window-size "$WINDOW_SIZE" --fst-window-step "$WINDOW_STEP" --out "$out"
        [[ -s "${out}.windowed.weir.fst" ]] || { echo "ERROR: missing output for $n1 vs $n2" >&2; exit 3; }
    done
done
