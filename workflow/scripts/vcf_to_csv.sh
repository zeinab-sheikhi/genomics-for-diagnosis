#!/bin/bash
# vcf_to_csv.sh

VCF_FILE="$1"
OUTPUT_CSV="$2"
LOG_FILE="$3"

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOG_FILE"
}

log_message "Starting VCF to CSV conversion"
log_message "Input VCF: $VCF_FILE"
log_message "Output CSV: $OUTPUT_CSV"

mkdir -p "$(dirname "$OUTPUT_CSV")"

if [[ ! -f "$VCF_FILE" ]]; then
    log_message "ERROR: VCF file not found: $VCF_FILE"
    exit 1
fi

TOTAL_VARIANTS=$(bcftools view -H "$VCF_FILE" | wc -l)
log_message "Total variants found: $TOTAL_VARIANTS"

if [[ $TOTAL_VARIANTS -eq 0 ]]; then
    log_message "No structural variants found"
    echo "CHROM,START,END,SIZE,QUAL,FILTER" > "$OUTPUT_CSV"
    exit 0
fi

log_message "Extracting variant information..."
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%QUAL\t%FILTER\t%INFO/SVTYPE\n' "$VCF_FILE" | \
awk -v log_file="$LOG_FILE" '
BEGIN {
    OFS = ","
    count = 0
}
{
    count++
    chrom = $1
    start = $2
    end_field = $3
    svlen_field = $4 
    qual_field = $5
    filter_field = $6
    svtype = $7
    
    # Calculate END position
    if (end_field != ".") {
        end = end_field
    } else if (svlen_field != ".") {
        svlen = svlen_field
        if (svlen < 0) svlen = -svlen  # Make positive
        end = start + svlen
    } else {
        end = start + 1  # Default for unknown
    }
    
    # Calculate SIZE
    if (svlen_field != ".") {
        size = svlen_field
        if (size < 0) size = -size  # Make positive
    } else {
        size = end - start
    }
    
    # Format QUAL
    if (qual_field != ".") {
        qual = sprintf("%.1f", qual_field)
    } else {
        qual = "0.0"
    }
    
    # Format FILTER
    if (filter_field == "." || filter_field == "") {
        filter = "PASS"
    } else {
        filter = filter_field
    }

    if (svtype == "." || svtype == "") {
        svtype  "UNK" 
    } else {
        svtype = svtype
    }
    
    # Output CSV row
    print chrom, start, end, size, qual, filter, svtype
}
END {
    print "Processed " count " variants" >> log_file
}' > temp_sv_data.csv

log_message "Sorting variants by genomic positions..."
(
    echo "CHROM,START,END,SIZE,QUAL,FILTER,SVTYPE"
    sort -t',' -k1,1V -k2,2n temp_sv_data.csv
) > "$OUTPUT_CSV"

rm -f temp_sv_data.csv

FINAL_COUNT=$(tail -n +2 "$OUTPUT_CSV" | wc -l)
log_message "Successfully exported $FINAL_COUNT variants to $OUTPUT_CSV"
log_message "VCF to CSV conversion completed successfully"

echo "Conversion completed: $FINAL_COUNT variants exported"