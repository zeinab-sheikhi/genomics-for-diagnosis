#!/bin/bash
# annotate.sh
# Map structural variants to genes using overlap analysis

SV_BED="$1"
GENES_BED="$2"
OUTPUT_FILE="$3"

if [ $# -ne 3 ]; then
    echo "Usage: $0 <sv_bed> <genes_bed> <output_file>"
    exit 1
fi

if [ ! -f "$SV_BED" ]; then
    echo "Error: SV BED file not found: $SV_BED"
    exit 1
fi

if [ ! -f "$GENES_BED" ]; then
    echo "Error: Genes BED file not found: $GENES_BED"
    exit 1
fi

if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools is required but not found in PATH"
    echo "Please install bedtools in your conda environment"
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_FILE")"
TMP_FILE=$(mktemp)

echo "Mapping structural variants to genes..."

bedtools intersect -a "$SV_BED" -b "$GENES_BED" -wa -wb > "$TMP_FILE"

if [ $? -ne 0 ]; then
    echo "Error: bedtools intersect failed"
    rm -f "$TMP_FILE"
    exit 1
fi


sort -k1,1 -k2,2n -k3,3n "$TMP_FILE" \
| bedtools groupby -g 1-6 -c 10 -o distinct -delim ';' \
> "$OUTPUT_FILE"

rm -f "$TMP_FILE"

# Report results
overlap_count=$(wc -l < "$OUTPUT_FILE")
total_svs=$(wc -l < "$SV_BED")

echo "Results:"
echo "- Total SVs: $total_svs"
echo "- SVs overlapping genes: $overlap_count"
echo "- Output saved to: $OUTPUT_FILE"