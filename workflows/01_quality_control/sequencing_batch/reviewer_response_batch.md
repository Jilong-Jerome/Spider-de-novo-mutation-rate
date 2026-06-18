# Response to reviewer: sample pooling / processing in sequencing batches

*Three paste-ready blocks: (a) point-by-point reply, (b) Methods/Supplementary
paragraph, (c) figure/table captions. The sequencing provider is left unnamed;
the text describes the platform, flowcells and lanes only. Replace "Table Sx /
Figure Sx" with the final numbering.*

---

## (a) Point-by-point reply

> *"…missing details regarding exactly how samples were pooled / processed
> through machine(s) and in groups when sequencing was performed… determine
> if/how samples were processed in batches, and what samples were processed
> together… confirm there are not acute differences between the 'pairs' on the
> phylogenetic tree, and between parents vs. probands…"*

We thank the reviewer and have added the complete batching inventory to the
supplement (**Table Sx**, per-individual; **Figure Sx**, a coverage heatmap of
every individual × flowcell-lane) together with the batch-structure analysis
below.

**One platform.** All 271 individuals (7 species; 65 parents and 206 probands
across 34 families) were whole-genome sequenced on a **single platform,
DNBSEQ-G400**. The differing flowcell-ID prefixes that appear in the inventory
(V300…, V350…, DP8400…, FP100…) are **flowcell/run identifiers on this one
platform** — not different platforms, instruments or library-preparation kits.
There is therefore no platform or kit difference between any samples, and in
particular none between the phylogenetic pairs or between parents and probands.

**The only batch axis is the flowcell, and flowcell batches do not line up with
sociality.** At the per-base depth used here, sequencing was distributed
sample- and family-wise across **73 flowcells, each carrying only a few samples**
(median 4 individuals and ~2 families per flowcell; range 1–17 individuals),
with each family itself spread across a median of 4 flowcells (range 2–8). In
other words, samples were *not* pooled into a small number of large batches that
could coincide with the social/subsocial contrast. Indeed, **62 of the 73
flowcells functionally contain only a single species and 33 contain only a single
family** — i.e. the great majority of flowcell batches are nested within one
biological unit and so cannot, even in principle, separate social from subsocial
samples — while of the 11 flowcells that do mix species, 8 mix *across* the
social/subsocial divide and the other 3 combine two social species (*S. dumicola*
and *S. mimosarum*); not a single flowcell isolates social from subsocial.
Concretely:

1. **8 of the 73 flowcells carried both social and subsocial individuals** — the
   same physical flowcell contained both classes. For example flowcell
   `DP8400011737BL` carried *S. dumicola* and *S. mimosarum* (social) together
   with *S. tentoriicola* (subsocial); flowcells `DP8400011312BL`,
   `DP8400011454TL` and `DP8400011678BR` each carried *S. dumicola* (social) with
   *S. tentoriicola* (subsocial); and flowcells `V300088307`, `V300088321`,
   `V300088752` and `V300101732` each carried *S. lineatus* (subsocial) with
   *S. sarasinorum* (social). There is thus **no flowcell batch dedicated to
   "social" or to "subsocial" samples**.

2. **Flowcell assignment is unstructured with respect to parent vs. proband as
   well.** Of the 73 flowcells, 5 are parent-only, 34 are proband-only and **34
   carry both parents and probands**, and in **20 of the 34 families a parent and
   a proband of the same trio were co-sequenced on a shared flowcell** (the
   remainder are spread across flowcells because occupancy is sparse). There is
   thus no "parents-in-one-batch, probands-in-another" structure that could bias
   the within-trio de novo mutation comparison.

3. **Parents and probands are sequenced to comparable depth** (median ≈83×
   for both; mean ≈86× parents vs. ≈91× probands; Mann–Whitney U *p* ≈ 0.33, i.e.
   no significant difference, and if anything probands are marginally *better*
   covered, so offspring are not under-powered for de novo mutation detection).

Because de novo mutations are called *within each family* (a proband against its
own two parents) on identical chemistry, and because flowcell identity is
essentially haphazard with respect to both sociality and role, flowcell-level
technical variation cannot generate a systematic difference between the social
and subsocial lineages. We have added a summary of this batching scheme to the
Methods and the full inventory to the Supplementary Data (text below).

---

## (b) Methods / Supplementary paragraph

**Sequencing batches and processing.** All 271 individuals (65 parents and 206
probands across seven *Stegodyphus* species and 34 families) were whole-genome
sequenced on a single platform (DNBSEQ-G400). Each individual was sequenced on
one or more flowcells, with several individuals multiplexed per lane; the
flowcell-ID prefixes recorded in the inventory (V300…/V350…/DP8400…/FP100…) are
run/flowcell identifiers on this one platform and do not denote different
platforms or library kits. Samples were distributed across 73 flowcells, each
carrying only a few samples (median 4 individuals and ~2 families per flowcell),
and individual families were spread across a median of four flowcells; sequencing
was therefore not organised into large pooled batches. Most flowcells were
nested within a single biological unit (62 of 73 contained only one species and
33 only one family), and flowcell assignment was unstructured with respect to
sociality (8 flowcells carried both social and subsocial individuals) and to
parent/proband role (5 flowcells parent-only, 34 proband-only, 34 carrying both;
and in 20 of 34 families a parent and a proband of the same trio were co-sequenced
on a shared flowcell), so no flowcell batch coincides with the biological
comparisons. Two *S. mimosarum* males (`7.M` and `8.M`) are
each the sire of two families and so serve as a parent in two trios; the single
sequenced sample for each is used in both trios' within-family de novo mutation
calling (these two sires therefore appear once per family context in Figure Sx).
The complete per-individual mapping of species, family, parent/proband role,
flowcell IDs, lanes and total coverage is given in **Table Sx**, and shown
graphically as a per-individual × flowcell-lane coverage heatmap in **Figure
Sx**. Total coverage was comparable between parents and probands (median ≈83×
for both; Mann–Whitney *p* ≈ 0.33).

---

## (c) Supplementary captions

**Table Sx. Per-individual sequencing-batch inventory.** Each row is one
individual, giving its species, family, role (parent vs. proband), sociality, the
sequencing platform (DNBSEQ-G400 for all individuals), the specific flowcell
ID(s) and lane(s) on which it was sequenced, the number of flowcells/lane-level
records, and total genome-wide coverage. The accompanying summary tables report
the flowcell-batch structure (individuals and families per flowcell; flowcells
per family) and the flowcells that carried both social and subsocial individuals,
showing that no flowcell batch is aligned with sociality or with parent/proband
role. (Source: `batch_summary_table.tsv` / `batch_summary_table.md`;
batch-structure statistics in `batch_association_stats.txt`.)

**Figure Sx. Sequencing scheme as a coverage heatmap.** Rows are individuals
(grouped by sociality, species and family); columns are unique flowcell × lane
combinations on the single sequencing platform (DNBSEQ-G400). Cell colour encodes
approximate per-lane coverage; grey cells mark flowcell-lanes in which an
individual was not sequenced. The three left-hand annotation strips indicate,
from left to right, **Role** (parent vs. proband), **Sociality** (social vs.
subsocial) and **Species**. The sparse, scattered occupancy pattern — and the
absence of any block of columns dedicated to social vs. subsocial individuals —
shows that samples were not pooled into sociality-aligned sequencing batches.
(The figure has 273 rows rather than 271: two *S. mimosarum* males are each the
shared sire of two families and are shown once per family context, though each is
a single sequenced sample.) (Source: `sequencing_heatmap.pdf`.)

---

## Generated artifacts in this directory

```text
batch_summary_table.tsv       per-individual table (supplement)
batch_summary_table.md        table + flowcell-batch summaries + cross-sociality flowcells
sequencing_heatmap.pdf        Figure Sx (individual x flowcell-lane, with Role strip)
batch_association_stats.txt   platform note, flowcell-batch structure, cross-sociality
                              flowcells, coverage Mann-Whitney U
```
