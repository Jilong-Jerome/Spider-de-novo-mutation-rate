# Kinship Identification Workflow

Infer the **pairwise kinship relationship** between every pair of sequenced
individuals from genome-wide GATK-called genotypes, and in particular
distinguish **full-sibling** from **parent-offspring** relationships.

Full-sibs and parent-offspring are both 1st-degree (kinship coefficient ≈ 0.25),
so kinship alone cannot tell them apart. The discriminator is **IBS0** (and the
PLINK IBD probabilities Z0/Z1/Z2): parent-offspring share one allele at every
locus → IBS0 ≈ 0, Z1 ≈ 1; full-sibs → IBS0 > 0, Z0 ≈ 0.25, Z1 ≈ 0.5, Z2 ≈ 0.25.

This project modernizes the earlier ad-hoc `chapter3_dnm/script/kinship_check/`
pipeline into the config-driven gwf pattern, adds **KING** (gold-standard for
PO vs FS), and validates inferred relationships against the ground truth encoded
in the sample names.

## Pipeline (per species)

```
{SP}_filter_snps  -> {SP}_prepare_plink  -> {SP}_plink_ibd   ┐
                                          -> {SP}_king_kinship ┴-> {SP}_classify -> {SP}_plot

{SP}_trio_filter  -> {SP}_trio_check -> {SP}_trio_plot          (parallel branch)
```

The KING/PLINK branch estimates relatedness for **all pairs**. The `trio_*` branch is a
separate, allele-frequency-free validation of the labeled parent-offspring trios (see
"Per-trio Mendelian validation" below) — added because the frequency-based estimators
break down under social-spider inbreeding (true PO pairs get nonsensical, even negative,
KING kinship). The trio branch is the authoritative answer for PO confirmation; KING/PLINK
is kept as a cross-check.

**`run_king_plink` toggle** (config, default `true`): when `false` the 6 KING/PLINK
targets are not registered and only the trio branch runs. MIM/BIC keep it `true` (already
completed); the other species (AFR, DUM, LIN, SAR, TEN) run **trio-only** (`false`).

1. **filter_snps** (`vcftools`, `python_phylo`) — biallelic, PASS-only, well-genotyped
   (`--minGQ`, `--max-missing`), **autosomal** SNPs with MAF > `maf_threshold`. →
   `steps/{SP}/{SP}_autosome_snps_maf.vcf`
2. **prepare_plink** (`python_phylo`, `plink`) — fill variant IDs, LD-prune
   (`indep-pairwise {ld_prune}`), `--make-bed`, then remap scaffold names to
   integer chromosome codes for KING. → `steps/{SP}/{SP}_pruned.{bed,bim,fam}` +
   `{SP}_pruned.chrommap.tsv`
3. **plink_ibd** (`plink`) — `--genome` IBD (Z0/Z1/Z2/PI_HAT), independent
   cross-check. → `steps/{SP}/{SP}_plink.genome`
4. **king_kinship** (`king_binary`) — `--kinship --ibs` (Kinship + IBS0) and
   `--related --degree 2` (InfType). → `steps/{SP}/{SP}_king.kin0`,
   `{SP}_king_rel.kin0`
5. **classify** (`python_phylo`) — merge KING + PLINK, threshold-classify, derive
   known truth, compare. → `steps/{SP}/{SP}_kinship_classified.tsv`
6. **plot** (`python_phylo`) — kinship-vs-IBS0 + Z0-vs-Z1 + confusion matrix. →
   `steps/{SP}/{SP}_kinship_validation.pdf`

### Per-trio Mendelian validation (trio branch)

7. **trio_filter** (`vcftools`, `python_phylo`) — independent of `filter_snps`; same
   biallelic/PASS/`--minGQ` filter but `--max-missing {trio_max_missing}` (0.0 = keep all
   sites) and **no MAF cut** (`vcf_filter_maf_snps.py … -1`). → `{SP}_autosome_snps_nomaf.vcf`
8. **trio_check** (`python_phylo`, `trio_mendel_check.py`) — for each labeled `(F, M, S_k)`
   trio, scan sites where the two parents are **opposite homozygotes** (`0/0` vs `1/1`);
   a genuine offspring must be heterozygous there. Reports per offspring the het fraction
   (`consistency`) and a `verdict` (`consistent` / `INCONSISTENT` / `insufficient-sites`).
   → `{SP}_trio_mendel.tsv`
9. **trio_plot** (`python_phylo`, `plot_trio_mendel.py`) — consistency by family, power vs
   consistency, and a consistency histogram. → `{SP}_trio_mendel.pdf`

**combined_trio_plot** (`python_phylo`, `plot_combined_trio_mendel.py`) — a single
cross-species supplemental figure registered once (not per species) by
`combined_trio_plot()` in `workflow_sources.py`, downstream of every `{SP}_trio_check`.
One row per species (italic Latin binomials, alphabetical), two columns — consistency by
family + power-vs-consistency; the per-species histogram panel is dropped. Excludes
`BIC:family6` (the `--exclude SP:family` flag is repeatable). →
`steps/combined/all_species_trio_mendel.pdf`

Why: KING robust kinship is meaningless under this population's inbreeding (true PO pairs
came out at kinship −0.49…−13). The opposite-homozygote test uses only within-trio
Mendelian segregation — no allele frequencies — so it is unaffected. It confirms *both*
assigned parents are the true parents of an offspring; it does not by itself separate
full-sib from offspring (not needed here). Thresholds: `trio_consistency_min` (default
0.9), `trio_min_sites` (default 100), `trio_min_gq` (per-genotype GQ floor), `trio_min_dp`
(per-genotype depth floor, default 26 = the DNM callability depth). The depth floor guards
against low-depth allelic dropout miscalling a true het as homozygous — which would
otherwise fabricate a false opposite-homozygote parent site or a false offspring Mendelian
error. DP is read per-genotype inside `trio_mendel_check.py` (the no-MAF VCF carries
`GT:AD:DP:GQ:PL`), so changing it reruns only `trio_check`/`trio_plot`, not the filter pass.

## Autosomes only

Kinship/IBS0 must be computed on **autosomes only** — sex chromosomes distort
relatedness because males are hemizygous. `vcf_filter_maf_snps.py` keeps only
contigs named `<name>_<digits>`, so MIM's sex chromosomes `mim_X1`/`mim_X2` are
dropped (autosomes are `mim_1`…`mim_14`). `remap_bim_chrom.py` additionally
asserts no sex contig leaked into the KING marker set.

## Classification rules

Degree from KING robust kinship coefficient (cut-offs in config):

| kinship | class |
|---|---|
| > `kinship_dup` (0.354) | duplicate/MZ |
| `kinship_first`–`kinship_dup` (0.177–0.354) | 1st degree → split below |
| `kinship_second`–`kinship_first` | 2nd-degree |
| `kinship_third`–`kinship_second` | 3rd-degree |
| ≤ `kinship_third` (0.0442) | unrelated |

1st-degree split: **parent-offspring** if IBS0 ≤ `ibs0_po_max` (0.005), else
**full-sibling**. PLINK Z0 reported as cross-check.

## Ground truth from sample names

Sample IDs are `{SP}_family{N}_{role}_{...}` (`F_female`, `M_male`,
`S{k}_offspring`):

- different family → **unrelated**
- two parents (F & M) of same family → **unrelated** (mates/founders)
- parent vs offspring, same family → **parent-offspring**
- offspring vs offspring, same family → **full-sibling**

The classified TSV's `agree` column flags inferred-vs-known mismatches.

## Layout

```
workflow.py                       # entry point (glob configs)
configurations/mim_kinship.config.yaml
workflow_sources/
  workflow_sources.py             # builder (loops species) + combined_trio_plot (cross-species)
  workflow_templates.py           # 10 AnonymousTarget templates (6 KING/PLINK + 3 trio + 1 combined)
  vcf_filter_maf_snps.py          # reused: biallelic autosomal SNP filter (no-MAF mode via 3rd arg -1)
  vcf_fill_id.py                  # reused: fill "." variant IDs with chrom:pos
  remap_bim_chrom.py              # scaffold -> integer chrom for KING (asserts no sex chr)
  classify_kinship.py             # merge + threshold classify + truth comparison
  plot_kinship.py                 # validation figure
  trio_mendel_check.py            # per-trio opposite-homozygote Mendelian test
  plot_trio_mendel.py             # trio consistency figure (per species)
  plot_combined_trio_mendel.py    # combined all-species supplemental figure
steps/{SP}/                       # outputs
steps/combined/                   # cross-species combined figure
logs/{SP}_{step}.DONE             # sentinels
```

## Conda environments / tools

`vcftools`, `plink` (v1.9), `python_phylo` (pysam + pandas/matplotlib/numpy — serves
every Python step; there is no standalone `pysam` env), `gwf_new` (to drive gwf).
KING 2.2.7 is a static binary at `/home/jilong/software/king` (no conda env).

## Adding more species

Add entries under `species:` in the config (`{SP}: {vcf: <path>}`); the builder
registers the full 9-step pipeline (6 KING/PLINK + 3 trio) per species. GATK VCFs for all 7 species
(AFR, BIC, DUM, LIN, MIM, SAR, TEN) live at
`chapter3_dnm/steps/gatk/genotype/{SP}/{SP}_GATK.vcf.gz`.

## Running

```bash
cd kinship_identification
conda run -n gwf_new gwf status   # check the DAG
conda run -n gwf_new gwf run      # submit to SLURM
```
