#!/usr/bin/env python3
# bed2csv.py

import pandas as pd
import sys


def convert_bed_to_csv(annotated_bed, output_csv, large_sv_threshold=100000):
    """Convert annotated BED to CSV with gene lists and large SV flagging."""
    
    try:
        # Read annotated BED file
        df = pd.read_csv(annotated_bed, sep='\t', header=None,
                        names=['sv_chrom', 'sv_start', 'sv_end', 'sv_id', 'sv_size', 'sv_strand',
                               'gene_chrom', 'gene_start', 'gene_end', 'gene_name', 'gene_score', 'gene_strand'])
        
        # Parse SV information
        df['sv_start_1based'] = df['sv_start'] + 1  # Convert to 1-based
        
        # Extract SV details from ID
        sv_parts = df['sv_id'].str.split('_', expand=True)
        df['svtype'] = sv_parts[0]
        df['qual'] = sv_parts[2] if sv_parts.shape[1] > 2 else '.'
        df['filter'] = sv_parts[3] if sv_parts.shape[1] > 3 else '.'
        
        # Flag large SVs
        df['large_sv_flag'] = df['sv_size'].apply(lambda x: 'YES' if x >= large_sv_threshold else 'NO')
        
        # Group genes by SV
        result = df.groupby(['sv_chrom', 'sv_start_1based', 'sv_end', 'svtype', 'sv_size', 'qual', 'filter', 'large_sv_flag']).agg({
            'gene_name': lambda x: ';'.join(x.dropna().unique()) if x.dropna().size > 0 else 'INTERGENIC'
        }).reset_index()
        
        # Rename columns
        result.columns = ['CHROM', 'START', 'END', 'SVTYPE', 'SIZE', 'QUAL', 'FILTER', 'LARGE_SV_FLAG', 'OVERLAPPING_GENES']
        
        # Reorder columns
        result = result[['CHROM', 'START', 'END', 'SVTYPE', 'SIZE', 'QUAL', 'FILTER', 'OVERLAPPING_GENES', 'LARGE_SV_FLAG']]
        
        # Sort by chromosome and position
        result = result.sort_values(['CHROM', 'START'])
        
        # Save to CSV
        result.to_csv(output_csv, index=False)
        
        # Print statistics
        total_svs = len(result)
        large_svs = len(result[result['LARGE_SV_FLAG'] == 'YES'])
        gene_overlapping = len(result[result['OVERLAPPING_GENES'] != 'INTERGENIC'])
        
        print(f"Results:")
        print(f"- Total SVs: {total_svs}")
        print(f"- SVs overlapping genes: {gene_overlapping}")
        print(f"- Large SVs (≥{large_sv_threshold // 1000}kb): {large_svs}")
        print(f"- Output saved to: {output_csv}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python bed2csv.py <annotated_bed> <output_csv>")
        sys.exit(1)
    
    convert_bed_to_csv(sys.argv[1], sys.argv[2])
