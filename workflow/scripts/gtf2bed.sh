#!/bin/bash
# gtf_to_bed_chr21.sh
# Convert GENCODE GTF file to BED format for chromosome 21

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <gtf_file> [output_file]"
    echo "Example: $0 gencode.v48.basic.annotation.gtf.gz.Z chr21_genes.bed"
    exit 1
fi

GTF_FILE="$1"
OUTPUT_FILE="${2:-chr21_genes.bed}"

# Check if input file exists
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: File $GTF_FILE not found"
    exit 1
fi

echo "Extracting chromosome 21 genes from $GTF_FILE..."

# Convert GTF to BED
zcat "$GTF_FILE" | \
awk '$1 == "chr21" && $3 == "gene"' | \
while IFS=$'\t' read -r chr source feature start end score strand frame attributes; do
    # Extract gene_name only
    gene_name=$(echo "$attributes" | sed -n 's/.*gene_name "\([^"]*\)".*/\1/p')
    
    # Skip this gene if no gene_name found (should be rare)
    if [ -z "$gene_name" ]; then
        echo "Warning: No gene_name found for line: $chr:$start-$end" >&2
        continue
    fi

    echo -e "$chr\t$((start-1))\t$end\t$gene_name\t0\t$strand"
done | \
sort -k1,1 -k2,2n > "$OUTPUT_FILE"

gene_count=$(wc -l < "$OUTPUT_FILE")
echo "Extracted $gene_count genes to $OUTPUT_FILE"
echo "Done!"
