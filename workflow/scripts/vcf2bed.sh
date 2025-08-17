#!/bin/bash
# vcf_to_bed.sh
# Convert VCF file to BED format for bedtools analysis

VCF_FILE="$1"
OUTPUT_BED="$2"

# Check if input file is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <vcf_file> <output_bed>"
    echo "Example: $0 variants.vcf.gz variants.bed"
    exit 1
fi

# Check if input file exists
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: VCF file not found: $VCF_FILE"
    exit 1
fi

# Create output directory
mkdir -p "$(dirname "$OUTPUT_BED")"

echo "Converting VCF to BED format..."

# Extract SV information and convert to BED
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/SVLEN\t%QUAL\t%FILTER\n' "$VCF_FILE" | \
awk 'BEGIN{OFS="\t"} {
    chrom = $1
    pos = $2
    end_field = $3
    svtype = $4
    svlen_field = $5
    qual = $6
    filter = $7
    
    # Calculate END position
    if (end_field == ".") {
        if (svlen_field != ".") {
            svlen = (svlen_field < 0 ? -svlen_field : svlen_field)
            end = pos + svlen
        } else {
            end = pos + 1
        }
    } else {
        end = end_field
    }
    
    # Calculate SIZE
    if (svlen_field != ".") {
        size = (svlen_field < 0 ? -svlen_field : svlen_field)
    } else {
        size = end - pos
    }
    
    # Handle missing values
    if (qual == ".") qual = "0"
    if (filter == ".") filter = "PASS"
    if (svtype == ".") svtype = "UNK"
    
    # Create SV identifier
    sv_id = svtype "_" size "_" qual "_" filter
    
    # BED format: chrom, start(0-based), end, name, score, strand
    print chrom, (pos-1), end, sv_id, size, "+"
}' > "$OUTPUT_BED"

sv_count=$(wc -l < "$OUTPUT_BED")
echo "Converted $sv_count structural variants to BED format"
echo "Output saved to: $OUTPUT_BED"
