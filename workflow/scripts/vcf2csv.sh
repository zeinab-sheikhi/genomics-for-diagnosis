#!/bin/bash
# vcf2csv.sh

VCF_FILE="$1"

OUTPUT_CSV="$2"

mkdir -p "$(dirname "$OUTPUT_CSV")"

if [[ ! -f "$VCF_FILE" ]]; then
    echo "ERROR: VCF file not found: $VCF_FILE" >&2
    exit 1
fi

TOTAL_VARIANTS=$(bcftools view -H "$VCF_FILE" | wc -l)
echo "Total variants found: $TOTAL_VARIANTS"

if [[ $TOTAL_VARIANTS -eq 0 ]]; then
    echo "No structural variants found" >&2
    echo "CHROM,START,END,SIZE,QUAL,FILTER,SVTYPE" > "$OUTPUT_CSV"
    exit 0
fi

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%QUAL\t%FILTER\t%INFO/SVTYPE\n' "$VCF_FILE" | \
awk '
BEGIN { ="," }
{
    chrom=$1
    start=$2
    end_field=$3
    svlen_field=$4 
    qual_field=$5
    filter_field=$6
    svtype=$7
    
    # Calculate END
    if (end_field != ".") {
        end = end_field
    } else if (svlen_field != ".") {
        svlen = (svlen_field < 0 ? -svlen_field : svlen_field)
        end = start + svlen
    } else {
        end = start + 1
    }
    
    # SIZE
    if (svlen_field != ".") {
        size = (svlen_field < 0 ? -svlen_field : svlen_field)
    } else {
        size = end - start
    }
    
    # QUAL
    if (qual_field != ".") {
        qual = sprintf("%.1f", qual_field)
    } else {
        qual = "0.0"
    }
    
    # FILTER
    filter = (filter_field == "." || filter_field == "" ? "PASS" : filter_field)
    
    # SVTYPE
    if (svtype == "." || svtype == "") {
        svtype = "UNK"
    }
    
    print chrom, start, end, size, qual, filter, svtype
}' | sort -t',' -k1,1V -k2,2n | \
(
    echo "CHROM,START,END,SIZE,QUAL,FILTER,SVTYPE"
    cat
) > "$OUTPUT_CSV"

FINAL_COUNT=$(tail -n +2 "$OUTPUT_CSV" | wc -l)
echo "Exported $FINAL_COUNT variants to $OUTPUT_CSV"
