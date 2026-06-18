#!/usr/bin/env python3
"""
Compare expression of testis-specific genes between samples.

Analyzes TPM values to determine if candidate samples show
testis-specific expression patterns.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def load_testis_genes(gene_file: Path) -> set:
    """Load MIM testis homolog gene IDs."""
    genes = set()
    with open(gene_file) as f:
        for line in f:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return genes


def analyze_expression(tpm_matrix: pd.DataFrame, testis_genes: set,
                       control_samples: list, candidate_samples: list,
                       tpm_threshold: float = 1.0) -> dict:
    """Analyze expression patterns of testis-specific genes.

    Returns analysis results dictionary.
    """
    results = {}

    # Identify testis genes present in TPM matrix
    all_genes = set(tpm_matrix.index)
    testis_genes_present = testis_genes & all_genes
    background_genes = all_genes - testis_genes_present

    results['total_genes'] = len(all_genes)
    results['testis_genes_in_matrix'] = len(testis_genes_present)
    results['background_genes'] = len(background_genes)

    # Calculate expression for each sample
    sample_results = {}
    all_samples = control_samples + candidate_samples

    for sample in all_samples:
        if sample not in tpm_matrix.columns:
            print(f"Warning: Sample {sample} not found in TPM matrix")
            continue

        tpm_vals = tpm_matrix[sample]

        # Testis genes expression
        testis_tpm = tpm_vals[tpm_vals.index.isin(testis_genes_present)]
        testis_expressed = (testis_tpm > tpm_threshold).sum()
        testis_total = len(testis_tpm)
        testis_mean_tpm = testis_tpm.mean()
        testis_median_tpm = testis_tpm.median()

        # Background genes expression
        bg_tpm = tpm_vals[tpm_vals.index.isin(background_genes)]
        bg_expressed = (bg_tpm > tpm_threshold).sum()
        bg_total = len(bg_tpm)
        bg_mean_tpm = bg_tpm.mean()

        # Proportion of expressed testis genes
        prop_testis = testis_expressed / testis_total if testis_total > 0 else 0
        prop_bg = bg_expressed / bg_total if bg_total > 0 else 0

        # Fold enrichment
        fold_enrichment = prop_testis / prop_bg if prop_bg > 0 else float('inf')

        # Fisher's exact test
        # 2x2 table: [testis_expressed, testis_not] vs [bg_expressed, bg_not]
        contingency = [
            [testis_expressed, testis_total - testis_expressed],
            [bg_expressed, bg_total - bg_expressed]
        ]
        odds_ratio, fisher_p = stats.fisher_exact(contingency, alternative='greater')

        sample_results[sample] = {
            'testis_expressed': testis_expressed,
            'testis_total': testis_total,
            'testis_prop': prop_testis,
            'testis_mean_tpm': testis_mean_tpm,
            'testis_median_tpm': testis_median_tpm,
            'background_expressed': bg_expressed,
            'background_total': bg_total,
            'background_prop': prop_bg,
            'background_mean_tpm': bg_mean_tpm,
            'fold_enrichment': fold_enrichment,
            'fisher_odds': odds_ratio,
            'fisher_pvalue': fisher_p,
            'sample_type': 'control' if sample in control_samples else 'candidate'
        }

    results['samples'] = sample_results
    return results


def generate_report(results: dict, output_file: Path):
    """Generate text report of analysis results."""
    with open(output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("TESTIS VERIFICATION REPORT\n")
        f.write("=" * 70 + "\n\n")

        f.write("SUMMARY\n")
        f.write("-" * 70 + "\n")
        f.write(f"Total genes in analysis: {results['total_genes']}\n")
        f.write(f"Testis-specific gene homologs: {results['testis_genes_in_matrix']}\n")
        f.write(f"Background genes: {results['background_genes']}\n\n")

        f.write("SAMPLE ANALYSIS\n")
        f.write("-" * 70 + "\n\n")

        for sample, data in results['samples'].items():
            sample_type = data['sample_type'].upper()
            f.write(f"Sample: {sample} [{sample_type}]\n")
            f.write(f"  Testis genes expressed (TPM > 1): {data['testis_expressed']} / {data['testis_total']} ({data['testis_prop']*100:.1f}%)\n")
            f.write(f"  Background genes expressed: {data['background_expressed']} / {data['background_total']} ({data['background_prop']*100:.1f}%)\n")
            f.write(f"  Mean TPM of testis genes: {data['testis_mean_tpm']:.2f}\n")
            f.write(f"  Median TPM of testis genes: {data['testis_median_tpm']:.2f}\n")
            f.write(f"  Fold enrichment: {data['fold_enrichment']:.2f}x\n")
            f.write(f"  Fisher's exact test p-value: {data['fisher_pvalue']:.2e}\n\n")

        # Interpretation
        f.write("INTERPRETATION\n")
        f.write("-" * 70 + "\n")

        # Get control and candidate stats
        control_props = [d['testis_prop'] for s, d in results['samples'].items() if d['sample_type'] == 'control']
        candidate_props = [d['testis_prop'] for s, d in results['samples'].items() if d['sample_type'] == 'candidate']

        if control_props and candidate_props:
            control_mean = np.mean(control_props)
            candidate_mean = np.mean(candidate_props)

            f.write(f"\nControl samples mean testis gene expression: {control_mean*100:.1f}%\n")
            f.write(f"Candidate samples mean testis gene expression: {candidate_mean*100:.1f}%\n")

            if candidate_mean > control_mean * 1.5:
                f.write("\nCONCLUSION: Candidate samples show ELEVATED expression of testis-specific\n")
                f.write("genes compared to ovary control, consistent with testis tissue.\n")
            elif candidate_mean > control_mean * 1.2:
                f.write("\nCONCLUSION: Candidate samples show MODERATELY ELEVATED expression of\n")
                f.write("testis-specific genes. Further investigation recommended.\n")
            else:
                f.write("\nCONCLUSION: Candidate samples do NOT show elevated expression of\n")
                f.write("testis-specific genes compared to control.\n")

        f.write("\n" + "=" * 70 + "\n")


def generate_heatmap(tpm_matrix: pd.DataFrame, testis_genes: set,
                     samples: list, output_file: Path):
    """Generate expression heatmap of testis-specific genes."""
    # Filter to testis genes
    testis_genes_present = [g for g in tpm_matrix.index if g in testis_genes]

    if not testis_genes_present:
        print("No testis genes found in matrix, skipping heatmap")
        return

    # Get TPM values for testis genes
    tpm_subset = tpm_matrix.loc[testis_genes_present, samples].copy()

    # Filter to expressed genes (at least one sample with TPM > 1)
    expressed_mask = (tpm_subset > 1).any(axis=1)
    tpm_expressed = tpm_subset[expressed_mask]

    if len(tpm_expressed) == 0:
        print("No expressed testis genes found, skipping heatmap")
        return

    # Log transform for visualization
    tpm_log = np.log2(tpm_expressed + 1)

    # Create heatmap
    fig, ax = plt.subplots(figsize=(10, max(8, len(tpm_log) * 0.15)))

    # Limit genes for readability
    if len(tpm_log) > 100:
        # Show top 100 most variable genes
        tpm_log = tpm_log.loc[tpm_log.var(axis=1).nlargest(100).index]

    sns.heatmap(tpm_log, cmap='viridis', ax=ax,
                xticklabels=True, yticklabels=len(tpm_log) <= 50,
                cbar_kws={'label': 'log2(TPM + 1)'})

    ax.set_xlabel('Sample')
    ax.set_ylabel('Testis-specific gene homolog')
    ax.set_title(f'Expression of Testis-Specific Gene Homologs\n({len(tpm_log)} genes shown)')

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Heatmap saved to {output_file}")


def generate_barplot(results: dict, output_file: Path):
    """Generate bar plot comparing testis gene expression across samples."""
    samples = list(results['samples'].keys())
    props = [results['samples'][s]['testis_prop'] * 100 for s in samples]
    colors = ['blue' if results['samples'][s]['sample_type'] == 'control' else 'red'
              for s in samples]

    fig, ax = plt.subplots(figsize=(8, 6))

    bars = ax.bar(samples, props, color=colors, edgecolor='black')

    ax.set_ylabel('% Testis Genes Expressed (TPM > 1)')
    ax.set_xlabel('Sample')
    ax.set_title('Expression of Testis-Specific Gene Homologs')

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', edgecolor='black', label='Control (Ovary)'),
                       Patch(facecolor='red', edgecolor='black', label='Candidate (Testis?)')]
    ax.legend(handles=legend_elements, loc='upper right')

    # Add value labels on bars
    for bar, prop in zip(bars, props):
        height = bar.get_height()
        ax.annotate(f'{prop:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')

    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"Bar plot saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Compare testis-specific gene expression between samples'
    )
    parser.add_argument('--tpm', required=True, type=Path,
                        help='TPM matrix file')
    parser.add_argument('--testis-genes', required=True, type=Path,
                        help='MIM testis homolog gene list')
    parser.add_argument('--control-samples', required=True, nargs='+',
                        help='Control sample names')
    parser.add_argument('--candidate-samples', required=True, nargs='+',
                        help='Candidate sample names')
    parser.add_argument('--output-dir', required=True, type=Path,
                        help='Output directory')
    parser.add_argument('--tpm-threshold', type=float, default=1.0,
                        help='TPM threshold for "expressed" (default: 1.0)')

    args = parser.parse_args()

    print(f"Loading TPM matrix from {args.tpm}...")
    tpm_matrix = pd.read_csv(args.tpm, sep='\t', index_col=0)
    print(f"  Loaded {len(tpm_matrix)} genes, {len(tpm_matrix.columns)} columns")

    print(f"Loading testis gene list from {args.testis_genes}...")
    testis_genes = load_testis_genes(args.testis_genes)
    print(f"  Loaded {len(testis_genes)} testis gene homologs")

    print("\nAnalyzing expression patterns...")
    results = analyze_expression(
        tpm_matrix, testis_genes,
        args.control_samples, args.candidate_samples,
        args.tpm_threshold
    )

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate report
    report_file = args.output_dir / 'testis_verification_report.txt'
    generate_report(results, report_file)
    print(f"\nReport written to {report_file}")

    # Generate visualizations
    all_samples = args.control_samples + args.candidate_samples
    heatmap_file = args.output_dir / 'expression_heatmap.png'
    generate_heatmap(tpm_matrix, testis_genes, all_samples, heatmap_file)

    barplot_file = args.output_dir / 'expression_barplot.png'
    generate_barplot(results, barplot_file)

    # Print summary to console
    print("\n" + "=" * 50)
    print("QUICK SUMMARY")
    print("=" * 50)
    for sample, data in results['samples'].items():
        pval_str = f"{data['fisher_pvalue']:.2e}"
        print(f"{sample}: {data['testis_prop']*100:.1f}% testis genes expressed "
              f"(fold={data['fold_enrichment']:.1f}x, p={pval_str})")

    return 0


if __name__ == '__main__':
    sys.exit(main())
