# DNM Spectrum Compare

SBS96 trinucleotide mutational spectrum analysis and per-site mutation rate estimation for germline and somatic de novo mutations across 7 spider species.

## Project Structure

```
DNM_spectrum_compare/
├── data/
│   ├── dnm/
│   │   ├── germline/{SP}_minDP_26_automate_classified.tsv      # Germline DNM calls (symlinks)
│   │   └── somatic/{SP}_final_mutation_table_no_clusters.tsv   # Somatic raw mutation tables (symlinks)
│   ├── callable_sites/
│   │   ├── germline/{SP}_callable_minDP26.tsv                  # Per (offspring, chrom) callable bp
│   │   └── somatic/{SP}_trinuc_counts.tsv                      # 32-pyrimidine-canonical trinuc opportunity totals
│   └── genome/{SP}_ncbi_chromosome.fa                          # Reference genomes (symlinks)
├── dnm_spectrum/                                                # gwf workflow: spectra + tests
│   ├── workflow.py                                              # Entry point
│   ├── configurations/spectrum.config.yaml                      # All parameters
│   └── workflow_sources/
│       ├── workflow_sources.py                                  # Workflow builder
│       ├── workflow_templates.py                                # GWF target templates
│       ├── extract_context.py                                   # Germline trinuc extraction + SBS96 assignment
│       ├── compute_spectrum.py                                  # Germline spectrum aggregation
│       ├── compute_somatic_spectrum.py                          # Local 96-cat somatic spectrum from raw table
│       ├── plot_spectrum.py                                     # SBS96 bar plot (6 panels, canonical colors)
│       ├── plot_comparison.py                                   # Germline vs somatic comparison + per-class test
│       └── ...                                                  # Test scripts (see Workflow Steps)
├── dnm_rate/                                                    # gwf workflow: per-site mutation rates
│   ├── workflow.py
│   ├── configurations/rate.config.yaml
│   └── workflow_sources/
│       ├── workflow_sources.py / workflow_templates.py
│       ├── callable_vcf_trinuc.py                               # Per-(offspring, chrom) trinuc count from callable VCF
│       ├── aggregate_callable_trinuc.py                         # Sum to species-level autosome trinuc table
│       ├── germline_rate_summary.py                             # Per-species germline rate (overall + 7 classes)
│       ├── somatic_rate_summary.py                              # Per-species somatic rate (overall + 7 classes)
│       ├── aggregate_rates.py                                   # Master TSV + per-trio + per-species somatic + fold-change
│       └── plot_rates.py                                        # Per-species/group rate plots + fold-change plot
├── steps/                                                       # Output
│   ├── genome_index/{SP}/
│   ├── annotated_dnm/{SP}/{SP}_dnm_annotated.tsv
│   ├── somatic_spectrum/{SP}/{SP}_somatic_spectrum.tsv
│   ├── spectrum/{group}/{group}_{chrom_group}_spectrum.tsv
│   ├── comparison/{group}/{group}_germline_vs_somatic.{pdf,test.tsv}
│   ├── tests/...
│   └── rate/
│       ├── germline_callable_trinuc/{SP}/...
│       ├── germline/{SP}/{SP}_germline_rate_summary.tsv
│       ├── somatic/{SP}/{SP}_somatic_rate_summary.tsv
│       ├── master/mutation_rate_master.tsv
│       ├── master/per_trio_germline_mutations_callable.tsv
│       ├── master/per_species_somatic_mutations_callable.tsv
│       ├── master/somatic_vs_germline_fold_change.tsv
│       └── plots/{per_species_rates,group_rates,fold_change}.pdf
└── logs/
```

## Species

AFR, BIC, DUM, LIN, MIM, SAR, TEN (7 spider species)

## Workflow Steps

### `dnm_spectrum/` — spectra and statistical tests
1. **`{SP}_index_genome`** — `samtools faidx` on genome FASTA (conda: `samtools117`)
2. **`{SP}_extract_context`** — Extract trinucleotide context from genome, apply strand collapsing to pyrimidine reference, assign SBS96 category (conda: `samtools117`)
3. **`{SP}_somatic_spectrum`** — Aggregate the raw somatic mutation table (already strand-collapsed via `oriented_ref/oriented_alt/oriented_context`) into a 96-category SBS spectrum locally; replaces the prior dependency on the external `somatic/` project (conda: `python_phylo`)
4. **`{group}_{chrom_group}_spectrum`** — Merge species per group, deduplicate shared DNMs, count 96 SBS categories (conda: `python_phylo`)
5. **`{group}_{chrom_group}_plot`** — Render SBS96 bar plot PDF with 6 mutation-type panels and canonical colors (conda: `python_phylo`)
6. **`{group}_comparison_plot`** — Render germline (autosome) vs merged-somatic SBS96 comparison PDF (proportions with bootstrap CI) plus a companion TSV of per-category two-sided Fisher's exact tests (BH-corrected at both 96- and 7-category levels); BH-significant categories are bolded on the plot x-axis (conda: `python_phylo`)
7. **`social_vs_subsocial_7category_test`** — Chi-square + per-category Fisher's exact test (BH-corrected) comparing social vs subsocial across 7 collapsed mutation categories, for germline and somatic; outputs TSV + annotated bar plot (conda: `python_phylo`)
8. **`social_vs_subsocial_deviation_interaction_test`** — Per-category Breslow-Day + global 3-way log-linear interaction test of whether germline drifts from somatic in different ways between social and subsocial species; outputs TSV + log-OR plot (conda: `python_phylo`)
9. **`pairwise_somatic_test`** — Per-pair (DUM/TEN, SAR/BIC, MIM/AFR) Fisher's exact + across-pair Cochran-Mantel-Haenszel test on somatic spectra to check whether per-category social-vs-subsocial differences are consistent in each closest social-subsocial pair; outputs TSV + per-category log-OR plot (conda: `python_phylo`)
10. **`baseline_somatic_test`** — Treat the merged subsocial somatic spectrum as the ancestral baseline (leave-one-out for subsocial species) and Fisher-test each species' somatic spectrum against it per category; heatmap of log2 OR plus social-species consistency bar (conda: `python_phylo`)
11. **`spectrum_phylo_mantel`** — Mantel permutation test of pairwise somatic-spectrum cosine distance (96- and 7-cat) vs pairwise phylogenetic (tree-topology) distance; non-significant r justifies the merged-subsocial baseline used in `baseline_somatic_test` (conda: `python_phylo`)

### `dnm_rate/` — per-site mutation rate estimation
1. **`{SP}_{offspring}_{chrom}_callable_trinuc`** — Stream the per-(offspring, chrom) callable VCF, fetch the trinucleotide context from the indexed genome at every callable site, canonicalize to pyrimidine middle, and write a 32-row trinuc-count TSV. One target per autosome VCF (parallelized, kept-offspring only) (conda: `samtools117`)
2. **`{SP}_germline_callable_trinuc`** — Sum the per-(offspring, chrom) 32-trinuc TSVs into a species-level autosome callable trinuc table (conda: `python_phylo`)
3. **`{SP}_germline_rate_summary`** — Count autosomal germline DNMs per 7-category class + overall after exclusions and dedup; derive class-conditioned denominators from the species trinuc table; emit per-species germline rate TSV with Poisson 95% CI (conda: `python_phylo`)
4. **`{SP}_somatic_rate_summary`** — Count somatic DNMs per 7-category class + overall from the raw oriented mutation table; derive class-conditioned denominators from the somatic trinuc_counts file; emit per-species somatic rate TSV with Poisson 95% CI (conda: `python_phylo`)
5. **`mutation_rate_aggregate`** — Combine all per-species summaries with the per-(offspring, chrom) trinuc tables to produce: master rate TSV (per species / per group / overall × 8 classes), per-trio per-chromosome germline mutations + callable TSV, per-species somatic mutations + callable TSV, somatic-vs-germline fold-change TSV (conda: `python_phylo`)
6. **`mutation_rate_plot`** — Render per-species, per-group, and fold-change rate plots (conda: `python_phylo`)

## Configuration

All parameters in `dnm_spectrum/configurations/spectrum.config.yaml` (also consumed by `dnm_rate/`):
- **`species.{SP}.dnm_file`** — Germline DNM TSV (`data/dnm/germline/...`)
- **`species.{SP}.somatic_dnm_file`** — Raw somatic mutation table (`data/dnm/somatic/...`)
- **`species.{SP}.germline_callable_file`** — Per-(offspring, chrom) callable bp TSV
- **`species.{SP}.somatic_trinuc_counts_file`** — Pyrimidine-canonical trinuc opportunity totals
- **`species.{SP}.x_chromosomes`** — X chromosome names; produces separate autosome and x_chromosome spectra; rate workflow computes autosomes only
- **`species.{SP}.exclude_trios`** — Offspring IDs to remove (BIC family2/3 for hyper-mutation; SAR family3 for low callability; both excluded for spectrum and rate calculations)
- **`somatic_spectra`** — Paths to locally produced 96-cat somatic spectra (output of `{SP}_somatic_spectrum`)
- **`germline_callable_vcf_root`** — Root of per-(species, offspring, chrom) callable VCFs (consumed by `dnm_rate`)
- **`spectrum_groups`** — Named groups of species to merge (e.g., social, subsocial, all_species)

## Running

```bash
# spectra + tests
cd dnm_spectrum
conda run -n gwf_new gwf status    # check targets
conda run -n gwf_new gwf run       # submit jobs

# mutation rates
cd ../dnm_rate
conda run -n gwf_new gwf status
conda run -n gwf_new gwf run
```

## Key Implementation Details

- **Strand collapsing**: Mutations reported with pyrimidine reference (C or T). Purine refs are reverse-complemented along with their trinucleotide context. Somatic data is pre-collapsed via `oriented_ref/oriented_alt/oriented_context`.
- **Deduplication**: Shared germline DNMs among siblings (same chrom/pos/ref/alt/parents) counted once at the spectrum/group-rate level. Per-trio summaries report the un-deduped per-offspring count.
- **Class-conditioned denominators (7-class)**: `C>A`/`C>G` use sum of NCN canonical trinucs (middle = C); `C>T_CpG` uses NCG; `C>T_nonCpG` = `C-middle total − CpG`; `T>A`/`T>C`/`T>G` use sum of NTN; `overall` uses all 32. Built from the somatic `trinuc_counts.tsv` for somatic; from the per-(offspring, chrom) callable VCF trinuc tables for germline (exact, not proxy).
- **Output**: 96-row TSV per group per chromosome group; per-species rate TSVs with Poisson 95% CI; master rate TSV with per-species, per-group, fold-change rows.
