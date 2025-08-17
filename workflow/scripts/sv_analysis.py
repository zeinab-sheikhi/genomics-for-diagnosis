#!/usr/bin/env python3
"""
Simple SV analysis script
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import sys


def parse_args():
    parser = argparse.ArgumentParser(description='Simple SV analysis')
    parser.add_argument('--vcf', required=True, help='Input VCF file')
    parser.add_argument('--csv', required=True, help='Input CSV file') 
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--outdir', required=True, help='Output directory')
    return parser.parse_args()


def load_data(csv_file):
    """Load SV data from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        print(f"Loaded {len(df)} structural variants")
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)


def create_genome_overview(df, output_file, sample_name):
    """Create basic genome overview plot."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle(f'SV Analysis Overview - {sample_name}', fontsize=16)
    
    # SV type distribution
    sv_counts = df['SVTYPE'].value_counts()
    ax1.pie(sv_counts.values, labels=sv_counts.index, autopct='%1.1f%%')
    ax1.set_title('SV Type Distribution')
    
    # Position distribution
    ax2.hist(df['START'], bins=50, alpha=0.7, color='skyblue')
    ax2.set_xlabel('Genomic Position')
    ax2.set_ylabel('Count')
    ax2.set_title('SV Position Distribution')
    
    # Size distribution
    ax3.hist(np.log10(df['SIZE'] + 1), bins=50, alpha=0.7, color='lightgreen')
    ax3.set_xlabel('Log10(Size + 1)')
    ax3.set_ylabel('Count')
    ax3.set_title('SV Size Distribution')
    
    # Quality distribution
    ax4.hist(df['QUAL'], bins=50, alpha=0.7, color='salmon')
    ax4.set_xlabel('Quality Score')
    ax4.set_ylabel('Count') 
    ax4.set_title('Quality Score Distribution')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Genome overview saved to: {output_file}")


def create_size_analysis(df, output_file, sample_name):
    """Create size distribution analysis."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f'SV Size Analysis - {sample_name}', fontsize=16)
    
    # Box plot by type
    sns.boxplot(data=df, x='SVTYPE', y='SIZE', ax=ax1)
    ax1.set_yscale('log')
    ax1.set_title('Size Distribution by SV Type')
    ax1.tick_params(axis='x', rotation=45)
    
    # Size categories
    size_cats = pd.cut(df['SIZE'], 
                      bins=[0, 100, 1000, 10000, 100000, float('inf')],
                      labels=['<100bp', '100bp-1kb', '1kb-10kb', '10kb-100kb', '>100kb'])
    size_counts = size_cats.value_counts()
    ax2.bar(range(len(size_counts)), size_counts.values, color='lightcoral')
    ax2.set_xticks(range(len(size_counts)))
    ax2.set_xticklabels(size_counts.index, rotation=45)
    ax2.set_ylabel('Count')
    ax2.set_title('SV Count by Size Category')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Size analysis saved to: {output_file}")


def create_quality_analysis(df, output_file, sample_name):
    """Create quality analysis plots."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f'SV Quality Analysis - {sample_name}', fontsize=16)
    
    # Quality by filter status
    pass_qual = df[df['FILTER'] == 'PASS']['QUAL']
    lowqual = df[df['FILTER'] != 'PASS']['QUAL']
    
    ax1.hist(pass_qual, bins=30, alpha=0.7, label=f'PASS (n={len(pass_qual)})', color='lightgreen')
    ax1.hist(lowqual, bins=30, alpha=0.7, label=f'LowQual (n={len(lowqual)})', color='lightcoral')
    ax1.set_xlabel('Quality Score')
    ax1.set_ylabel('Count')
    ax1.set_title('Quality by Filter Status')
    ax1.legend()
    
    # Quality vs Size scatter
    ax2.scatter(df['SIZE'], df['QUAL'], alpha=0.6, s=20)
    ax2.set_xlabel('SV Size (bp)')
    ax2.set_ylabel('Quality Score')
    ax2.set_title('Quality vs Size')
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Quality analysis saved to: {output_file}")


def main():
    args = parse_args()
    
    # Load data
    df = load_data(args.csv)
    
    # Create output directory
    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    genome_plot = f"{args.outdir}/{args.sample}_genome_overview.png"
    size_plot = f"{args.outdir}/{args.sample}_size_distribution.png"
    quality_plot = f"{args.outdir}/{args.sample}_quality_analysis.png"
    
    create_genome_overview(df, genome_plot, args.sample)
    create_size_analysis(df, size_plot, args.sample)
    create_quality_analysis(df, quality_plot, args.sample)


if __name__ == "__main__":
    main()
