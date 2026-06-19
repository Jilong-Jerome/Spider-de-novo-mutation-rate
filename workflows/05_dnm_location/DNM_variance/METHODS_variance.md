# Supplemental Methods: Individual-level test for a difference in de novo mutation-rate variability between social and subsocial species

## Overview

To ask whether the *variability* of the autosomal de novo mutation (DNM) rate
differs between social and subsocial *Stegodyphus*, we quantified, for each
species, how much the per-offspring mutation rate varies among the sequenced
trios, and then compared this variability between the two social systems against
a Poisson simulation null. Variability was summarised by the coefficient of
variation (CV = standard deviation / mean) of the per-individual rates, which
makes species with different mean rates directly comparable. The analysis was
applied to all seven species (three social: *S. dumicola* (DUM),
*S. mimosarum* (MIM), *S. sarasinorum* (SAR); four subsocial:
*S. tentoriicola* (TEN), *S. bicolor* (BIC), *S. africanus* (AFR),
*S. lineatus* (LIN)) using a configuration-driven Grid Workflow (`gwf`) pipeline
on a SLURM cluster. The analysis is reported at the **individual (trio) level**:
a parallel family-level test was not pursued, because with only three to four
families per species it has too little power to discriminate between competing
levels of inter-family variability. Species names, sex-chromosome labels,
excluded trios, the social groupings and the sister-species pairs are taken from
configuration files rather than hard-coded.

## De novo mutation set and family exclusions

The analysis was restricted to autosomes: the sex chromosomes listed per species
in the configuration (`x_chromosomes`, e.g. `afr_X1`, `afr_X2`; *S. lineatus* has
three) were removed before any computation. Within each offspring the DNM calls
were deduplicated to unique mutational events on the key
(chromosome, position, reference allele, alternate allele), and the unit of
analysis is the individual offspring (trio).

Four families were excluded for biological or technical reasons, configured as
`excluded_trios` (or removed upstream during data preparation):

- **Hypermutation.** *S. bicolor* families 2 and 3 carried anomalously elevated
  per-family rates (roughly three- to eight-fold the other families) and were
  removed; *S. bicolor* therefore contributes 18 rather than 30 individuals.
- **Zero calls / low callability.** *S. sarasinorum* family 3 (no detected DNMs
  over a low callable footprint) was removed, leaving 24 individuals; and
  *S. mimosarum* family 2, which was absent from the input DNM and callable-site
  data, leaving *S. mimosarum* with 18 individuals.

## Per-individual rate and per-species variation

For every offspring the per-base autosomal callable-site counts were summed to a
callable total (sex-chromosome records skipped), and the per-individual mutation
rate was computed as `unique DNMs / (2 × callable sites)`. The species mean rate was
estimated by pooling across all retained offspring as
`total unique DNMs / (2 × total callable sites)`; the factor of two accounts for
diploidy (each callable site contributes two alleles that could carry a new
mutation), giving a per-base-pair, per-generation rate with no further scaling. The
variability of a species was summarised by the coefficient of variation of its
per-individual rates, `CV = std(rates) / mean(rates)`, using the population
standard deviation (`ddof = 0`); the corresponding inter-individual variance was
also recorded.

## Poisson null model and test statistics

The null hypothesis is that, within a species, all offspring share that species'
single mean mutation rate and that their observed counts differ only through
Poisson sampling. For each simulation, every offspring's DNM count was drawn from
a Poisson distribution with mean `λ = 2 × (offspring callable sites) × (species mean rate)`,
divided by `2 × (offspring callable sites)` to give a simulated per-individual
rate, and the species CV was recomputed exactly as for the observed data. (The
factor of two cancels the diploid divisor in the mean rate, so the expected counts —
and hence the null distribution and all p-values — are identical to a per-callable-site
parameterisation; only the absolute rate scale changes.) This was
repeated for 10,000 simulations using NumPy's `default_rng` generator initialised
with a fixed seed (42), giving each species a null distribution of CV values that
embodies pure sampling noise around its own mean rate.

Two complementary social-versus-subsocial tests were then performed:

- **Global test.** The observed statistic is
  `mean CV(subsocial) − mean CV(social)`. Its null distribution is formed by taking
  the same difference across the per-species simulated CVs, and the one-tailed
  upper-tail empirical p-value is the fraction of simulations whose difference is
  greater than or equal to the observed difference (i.e. the evidence that
  subsocial species are *more* variable than social species).
- **Phylogenetically controlled sister-pair test.** Each subsocial species was
  paired with its closest social relative (bic:sar, afr:mim, ten:dum;
  *S. lineatus* has no social sister in this dataset and is omitted from the
  pairwise test). For each pair the observed `ΔCV = CV(subsocial) − CV(social)` was
  compared to the simulated null with a one-tailed upper p-value, and the number of
  pairs in which the subsocial member was the more variable was tabulated.

## Results

The seven species span a roughly three-fold range of individual-level CV, with no
systematic separation by social system; if anything the three social species sit
at the upper end.

| Species | Social system | N individuals | Mean rate | Individual CV |
|---|---|---|---|---|
| *S. africanus* (AFR)    | subsocial | 24 | 2.58e-09 | 0.456 |
| *S. bicolor* (BIC)      | subsocial | 18 | 3.81e-09 | 0.382 |
| *S. lineatus* (LIN)     | subsocial | 24 | 3.30e-09 | 0.292 |
| *S. tentoriicola* (TEN) | subsocial | 24 | 3.74e-09 | 0.611 |
| *S. dumicola* (DUM)     | social    | 24 | 2.18e-09 | 0.682 |
| *S. mimosarum* (MIM)    | social    | 18 | 1.48e-09 | 0.534 |
| *S. sarasinorum* (SAR)  | social    | 24 | 1.41e-09 | 0.815 |

The mean individual-level CV was 0.435 across the four subsocial species and 0.677
across the three social species. In the **global test** the observed difference
(subsocial − social) was −0.242, essentially identical to the null expectation
(null mean −0.232), giving a one-tailed p = 0.571: there is no evidence that
subsocial species are more variable, and the numerical trend is in the opposite
direction (social species slightly more variable).

The **sister-pair test** told the same story: in every one of the three
phylogenetically matched pairs the subsocial species was *less* variable than its
social relative, so zero of three pairs supported the subsocial-more-variable
hypothesis.

| Pair (subsocial:social) | CV subsocial | CV social | ΔCV | p (one-tailed upper) |
|---|---|---|---|---|
| bic:sar | 0.382 | 0.815 | −0.433 | 0.265 |
| afr:mim | 0.456 | 0.534 | −0.079 | 0.570 |
| ten:dum | 0.611 | 0.682 | −0.072 | 0.670 |

We therefore find **no significant difference in the individual-level variability
of the de novo mutation rate between social and subsocial *Stegodyphus***. The
observed differences are fully consistent with Poisson sampling noise, and what
little signal exists runs counter to the hypothesis that subsocial species harbour
more variable mutation rates.

## Reproducibility

The test was implemented as a `gwf` workflow on a SLURM cluster, with the
individual-level CV-comparison script
(`dnm_variance/workflow_sources/dnm_individual_cv_comparison.py`) and the final
output generator (`final_individual_variation/workflow_sources/make_final_individual_variation_outputs.py`)
pinned to a dedicated conda environment (`python_phylo`; Python 3 with NumPy,
pandas and PyYAML). Per-species inputs (DNM calls, callable-site files,
sex-chromosome labels and excluded trios) come from the species configuration
files, while the social groupings, sister-species pairs, number of simulations
(10,000) and random seed (42) come from `cv_comparison.config.yaml`. All inputs,
scripts, the per-species CV and variance tables, the per-simulation null draws and
the sister-pair test summaries are available in the project repository under
`steps/dnm_cv_comparison/` and `steps/final_individual_variation/`.
