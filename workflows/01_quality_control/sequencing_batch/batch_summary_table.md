# Sequencing batch summary

271 distinct sequenced individuals across 7 species and 34 families — 65 parents, 206 probands.

> Note: 2 sequenced males are each the shared sire of two families — `7.M` (S. mimosarum, families family2, family5); `8.M` (S. mimosarum, families family3, family4). They are one sample each (counted once above), but appear once per family context in Figure Sx, which therefore shows 273 rows.

**All individuals were sequenced on a single platform (DNBSEQ-G400).** `round_id` = sequencing flowcell/run; `lane` = lane within the flowcell. The differing flowcell-ID prefixes (V300…/V350…/DP8400…/FP100…) are run/flowcell identifiers on this one platform, not different platforms, instruments or library kits. The relevant "batch" unit is therefore the flowcell.

## Flowcell-batch structure

- Distinct flowcells: **73**
- Individuals per flowcell: median **4**, mean 5.3, range 1–17
- Families per flowcell: median **2**, mean 2.2, range 1–7
- Flowcells per family: median **4**, mean 4.7, range 2–8

  - individuals-per-flowcell distribution (count:flowcells): 1:7, 2:6, 3:5, 4:25, 5:4, 6:5, 7:6, 8:7, 9:1, 11:1, 12:2, 15:2, 16:1, 17:1
  - families-per-flowcell distribution: 1:33, 2:18, 3:10, 4:4, 5:5, 6:1, 7:2
  - flowcells-per-family distribution: 2:3, 3:6, 4:9, 5:4, 6:5, 7:6, 8:1

Sequencing was distributed sample-/family-wise across many small flowcells (not concentrated into a few large pooled batches), so no single batch could group all social trios apart from all subsocial trios.

## Flowcell composition (independence from the biological comparison)

How many of the 73 flowcells functionally contain samples from a single species / family / sociality class (nested within one biological unit) versus span more than one:

- Single-species flowcells: **62/73** (span >1 species: 11)
- Single-family flowcells: **33/73** (span >1 family: 40)
- Single-sociality flowcells: **65/73** (span both social & subsocial: 8)
- Single-individual flowcells: **7/73**

The two regimes both argue for independence of flowcell batch from sociality: most flowcells are nested within one species/family and so *cannot* separate social from subsocial samples, while of the 11 flowcells that do mix species, 8 mix *across* the social/subsocial divide (next section) and the other 3 combine species within one sociality class — so none isolates social from subsocial. In neither case does a flowcell batch line up with the biological comparison.

## Flowcells spanning both sociality classes

**8 of 73 flowcells carried BOTH social and subsocial individuals**, so flowcell batches do not segregate by sociality — the same flowcell physically contained both classes:

| Flowcell | Species on flowcell |
| --- | --- |
| DP8400011312BL | S. dumicola (social), S. tentoriicola (subsocial) |
| DP8400011454TL | S. dumicola (social), S. tentoriicola (subsocial) |
| DP8400011678BR | S. dumicola (social), S. tentoriicola (subsocial) |
| DP8400011737BL | S. dumicola (social), S. mimosarum (social), S. tentoriicola (subsocial) |
| V300088307 | S. lineatus (subsocial), S. sarasinorum (social) |
| V300088321 | S. lineatus (subsocial), S. sarasinorum (social) |
| V300088752 | S. lineatus (subsocial), S. sarasinorum (social) |
| V300101732 | S. lineatus (subsocial), S. sarasinorum (social) |

## Flowcell assignment is unstructured by role and by sociality

Because occupancy is sparse, the members of a family — including parents and probands — are themselves spread across several flowcells (median 4 flowcells per family), and the same flowcells are reused across families. Flowcell assignment therefore does not encode either biological variable:

- Flowcells carrying both parent and proband individuals: **34/73**
- Flowcells carrying both social and subsocial individuals: **8/73**

### Parent / proband flowcell composition

Breakdown of the 73 flowcells by the roles they functionally contain (parallels the species/family/sociality composition above):

- Parent-only flowcells: **5/73**
- Proband-only flowcells: **34/73**
- Mixed (both parents and probands): **34/73**
- Flowcells containing ≥1 parent: 39/73; ≥1 proband: 68/73
- Families in which a parent and a proband of that same family were co-sequenced on ≥1 shared flowcell: **20/34**

Parents and probands are therefore not segregated into separate sequencing batches: 34 of 73 flowcells carry both roles, and in 20 of 34 families a parent and a proband of the same trio were co-sequenced on one flowcell; the remainder are spread across flowcells. There is no "parents-in-one-batch, probands-in-another" structure that could bias the within-trio de novo mutation comparison.

Crucially, flowcell assignment is essentially haphazard with respect to the social/subsocial contrast (and to parent/proband role): there is no flowcell batch dedicated to one class, so flowcell-level technical variation cannot create a systematic social-vs-subsocial difference. De novo mutations are additionally called within each family on identical chemistry and filtered for run-specific artefacts.

## Per-individual table

| species | family | individual | role | sociality | platform | flowcell_ids | lanes | n_flowcells | n_lane_records | total_coverage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| S_dumicola | family1 | 51.S1 | proband | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 80.78 |
| S_dumicola | family1 | 52.S2 | proband | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 87.83 |
| S_dumicola | family1 | 53.S3 | proband | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 80.87 |
| S_dumicola | family1 | 54.S4 | proband | Social | DNBSEQ-G400 | DP8400011448TR | L01 | 1 | 1 | 101.42 |
| S_dumicola | family1 | 55.S5 | proband | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 77.8 |
| S_dumicola | family1 | 56.S6 | proband | Social | DNBSEQ-G400 | DP8400011448TR | L01 | 1 | 1 | 101.33 |
| S_dumicola | family1 | 41.F | parent | Social | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 97.72 |
| S_dumicola | family1 | 46.M | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 84.39 |
| S_dumicola | family2 | 57.S1 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 80.37 |
| S_dumicola | family2 | 58.S2 | proband | Social | DNBSEQ-G400 | DP8400011448TR | L01 | 1 | 1 | 112.38 |
| S_dumicola | family2 | 59.S3 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 77.64 |
| S_dumicola | family2 | 60.S4 | proband | Social | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 79.12 |
| S_dumicola | family2 | 61.S5 | proband | Social | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 68.8 |
| S_dumicola | family2 | 62.S6 | proband | Social | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 73.45 |
| S_dumicola | family2 | 42.F | parent | Social | DNBSEQ-G400 | DP8400011680BR | L01 | 1 | 1 | 137.41 |
| S_dumicola | family2 | 47.M | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 78.5 |
| S_dumicola | family3 | 63.S1 | proband | Social | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 78.85 |
| S_dumicola | family3 | 64.S2 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 91.29 |
| S_dumicola | family3 | 65.S3 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 80.29 |
| S_dumicola | family3 | 66.S4 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 79.85 |
| S_dumicola | family3 | 67.S5 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 81.64 |
| S_dumicola | family3 | 68.S6 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 81.44 |
| S_dumicola | family3 | 43.F | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 82.09 |
| S_dumicola | family3 | 48.M | parent | Social | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 100.37 |
| S_dumicola | family4 | 69.S1 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 83.17 |
| S_dumicola | family4 | 70.S2 | proband | Social | DNBSEQ-G400 | DP8400011448TR | L01 | 1 | 1 | 117.9 |
| S_dumicola | family4 | 71.S3 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 83.77 |
| S_dumicola | family4 | 72.S4 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 80.54 |
| S_dumicola | family4 | 73.S5 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 81.08 |
| S_dumicola | family4 | 74.S6 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 78.29 |
| S_dumicola | family4 | 44.F | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 68.33 |
| S_dumicola | family4 | 49.M | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 77.74 |
| S_dumicola | family5 | 75.S1 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 81.56 |
| S_dumicola | family5 | 76.S2 | proband | Social | DNBSEQ-G400 | DP8400011440TR;DP8400011680BR | L01 | 2 | 2 | 78.94 |
| S_dumicola | family5 | 77.S3 | proband | Social | DNBSEQ-G400 | DP8400011440TR | L01 | 1 | 1 | 63.5 |
| S_dumicola | family5 | 78.S4 | proband | Social | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 68.89 |
| S_dumicola | family5 | 79.S5 | proband | Social | DNBSEQ-G400 | DP8400011312BL;DP8400011737BL | L01 | 2 | 2 | 71.26 |
| S_dumicola | family5 | 80.S6 | proband | Social | DNBSEQ-G400 | DP8400011312BL;DP8400011737BL | L01 | 2 | 2 | 73.63 |
| S_dumicola | family5 | 45.F | parent | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 94.43 |
| S_dumicola | family5 | 50.M | parent | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 86.19 |
| S_mimosarum | family1 | 10.S2 | proband | Social | DNBSEQ-G400 | DP8400011456TL | L01 | 1 | 1 | 103.53 |
| S_mimosarum | family1 | 11.S3 | proband | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 99.69 |
| S_mimosarum | family1 | 12.S4 | proband | Social | DNBSEQ-G400 | DP8400011456TL;DP8400011737BL | L01 | 2 | 2 | 73.55 |
| S_mimosarum | family1 | 13.S5 | proband | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 70.52 |
| S_mimosarum | family1 | 14.S6 | proband | Social | DNBSEQ-G400 | DP8400011456TL | L01 | 1 | 1 | 93.87 |
| S_mimosarum | family1 | 9.S1 | proband | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 127.95 |
| S_mimosarum | family1 | 1.F | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 1 | 67.93 |
| S_mimosarum | family1 | 6.M | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 1 | 120.45 |
| S_mimosarum | family2 | 15.S1 | proband | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 125.17 |
| S_mimosarum | family2 | 16.S2 | proband | Social | DNBSEQ-G400 | DP8400011456TL | L01 | 1 | 1 | 102.86 |
| S_mimosarum | family2 | 17.S3 | proband | Social | DNBSEQ-G400 | DP8400011903TL | L01 | 1 | 1 | 153.68 |
| S_mimosarum | family2 | 18.S4 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 185.14 |
| S_mimosarum | family2 | 19.S5 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 190.53 |
| S_mimosarum | family2 | 20.S6 | proband | Social | DNBSEQ-G400 | DP8400011310BL | L01 | 1 | 1 | 73.22 |
| S_mimosarum | family2 | 2.F | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 1 | 137.38 |
| S_mimosarum | family2;family5 | 7.M | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 2 | 92.89 |
| S_mimosarum | family3 | 21.S1 | proband | Social | DNBSEQ-G400 | DP8400011903TL | L01 | 1 | 1 | 94.36 |
| S_mimosarum | family3 | 22.S2 | proband | Social | DNBSEQ-G400 | DP8400011310BL;DP8400011737BL | L01 | 2 | 2 | 72.39 |
| S_mimosarum | family3 | 23.S3 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 129.23 |
| S_mimosarum | family3 | 24.S4 | proband | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 117.09 |
| S_mimosarum | family3 | 25.S5 | proband | Social | DNBSEQ-G400 | DP8400011456TL | L01 | 1 | 1 | 120.57 |
| S_mimosarum | family3 | 26.S6 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 138.72 |
| S_mimosarum | family3 | 27.S7 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 131.7 |
| S_mimosarum | family3 | 28.S8 | proband | Social | DNBSEQ-G400 | DP8400011427BL | L01 | 1 | 1 | 176.16 |
| S_mimosarum | family3 | 3.F | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 1 | 87.49 |
| S_mimosarum | family3;family4 | 8.M | parent | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 2 | 91.94 |
| S_mimosarum | family4 | 29.S1 | proband | Social | DNBSEQ-G400 | DP8400011519TR | L01 | 1 | 1 | 97.58 |
| S_mimosarum | family4 | 30.S2 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 90.63 |
| S_mimosarum | family4 | 31.S3 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 89.09 |
| S_mimosarum | family4 | 32.S4 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 83.36 |
| S_mimosarum | family4 | 4.F | parent | Social | DNBSEQ-G400 | DP8400011455TL | L01 | 1 | 1 | 115.36 |
| S_mimosarum | family5 | 33.S1 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 78.35 |
| S_mimosarum | family5 | 34.S2 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 83.22 |
| S_mimosarum | family5 | 35.S3 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 85.22 |
| S_mimosarum | family5 | 36.S4 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 90.65 |
| S_mimosarum | family5 | 37.S5 | proband | Social | DNBSEQ-G400 | DP8400011680BR | L01 | 1 | 1 | 129.31 |
| S_mimosarum | family5 | 38.S6 | proband | Social | DNBSEQ-G400 | DP8400011903TL | L01 | 1 | 1 | 142.61 |
| S_mimosarum | family5 | 39.S7 | proband | Social | DNBSEQ-G400 | DP8400011680BR | L01 | 1 | 1 | 179.67 |
| S_mimosarum | family5 | 40.S8 | proband | Social | DNBSEQ-G400 | DP8400011452TL | L01 | 1 | 1 | 99.26 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S1 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 61.69 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S2 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 66.93 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S3 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 65.73 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S4 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 66.75 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S5 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 62.97 |
| S_sarasinorum | family1 | Ssar-A12-F-2-S6 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 93.49 |
| S_sarasinorum | family1 | Ssar-A12-F-2-F | parent | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 68.84 |
| S_sarasinorum | family1 | Ssar-A12-M | parent | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 81.25 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S1 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 77.44 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S2 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 76.19 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S3 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 76.2 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S4 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 75.57 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S5 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 76.09 |
| S_sarasinorum | family2 | Ssar-A17-F-4-S6 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 81.33 |
| S_sarasinorum | family2 | Ssar-A17-F-4-F | parent | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 75.44 |
| S_sarasinorum | family2 | Ssar-A17-M | parent | Social | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 95.27 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S1 | proband | Social | DNBSEQ-G400 | V300088299;V300088311;V300088769 | L1;L2;L3;L4 | 3 | 8 | 63.34 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S2 | proband | Social | DNBSEQ-G400 | V300088454;V300088899;V300101732 | L1;L2;L3;L4 | 3 | 12 | 66.03 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S3 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 68.77 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S4 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 61.32 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S5 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 71.08 |
| S_sarasinorum | family3 | Ssar-B11-F-3-S6 | proband | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 68.24 |
| S_sarasinorum | family3 | Ssar-B11-F-3-F | parent | Social | DNBSEQ-G400 | V300088454;V300088899 | L1;L2;L3;L4 | 2 | 8 | 64.07 |
| S_sarasinorum | family3 | Ssar-B11-M | parent | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 65.62 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S1 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 70.04 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S2 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 78.37 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S3 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 84.34 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S4 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 78.14 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S5 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 76.03 |
| S_sarasinorum | family4 | Ssar-B12-F-4-S6 | proband | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 81.72 |
| S_sarasinorum | family4 | Ssar-B12-F-4-F | parent | Social | DNBSEQ-G400 | V300088299;V300101733;V300102083 | L1;L2;L3;L4 | 3 | 8 | 79.75 |
| S_sarasinorum | family4 | Ssar-B12-M | parent | Social | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 93.09 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S1 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 63.32 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S2 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 66.76 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S3 | proband | Social | DNBSEQ-G400 | V300088905;V300088913 | L1;L2;L3;L4 | 2 | 8 | 71.73 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S4 | proband | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 74.4 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S5 | proband | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 70.16 |
| S_sarasinorum | family5 | Ssar-B25-F-5-S6 | proband | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 82.82 |
| S_sarasinorum | family5 | Ssar-B25-F-5-F | parent | Social | DNBSEQ-G400 | V300088751;V300088752 | L1;L2;L3;L4 | 2 | 6 | 73.49 |
| S_sarasinorum | family5 | Ssar-B25-M | parent | Social | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 94.97 |
| S_africanus | family1 | TB22_436_3_S1 | proband | Subsocial | DNBSEQ-G400 | V350231703 | L01;L02;L03;L04 | 1 | 4 | 83.1 |
| S_africanus | family1 | TB22_436_3_S2 | proband | Subsocial | DNBSEQ-G400 | V350231703 | L01;L02;L03;L04 | 1 | 4 | 83.17 |
| S_africanus | family1 | TB22_436_3_S3 | proband | Subsocial | DNBSEQ-G400 | V350231840 | L01;L02;L03;L04 | 1 | 4 | 83.13 |
| S_africanus | family1 | TB22_436_3_S4 | proband | Subsocial | DNBSEQ-G400 | V350231840 | L01;L02;L03;L04 | 1 | 4 | 83.17 |
| S_africanus | family1 | TB22_436_3_S5 | proband | Subsocial | DNBSEQ-G400 | V350231840 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family1 | TB22_436_3_S6 | proband | Subsocial | DNBSEQ-G400 | V350231840 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family1 | TB22_430-101_Ma | parent | Subsocial | DNBSEQ-G400 | V350231703 | L01;L02;L03;L04 | 1 | 4 | 83.13 |
| S_africanus | family1 | TB22_436_3_Fema | parent | Subsocial | DNBSEQ-G400 | V350231703 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family2 | TB22_440_5_S1 | proband | Subsocial | DNBSEQ-G400 | V350235689 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family2 | TB22_440_5_S2 | proband | Subsocial | DNBSEQ-G400 | V350235689 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family2 | TB22_440_5_S3 | proband | Subsocial | DNBSEQ-G400 | V350235693 | L01;L02;L03;L04 | 1 | 4 | 83.13 |
| S_africanus | family2 | TB22_440_5_S4 | proband | Subsocial | DNBSEQ-G400 | V350235693 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family2 | TB22_440_5_S5 | proband | Subsocial | DNBSEQ-G400 | V350235693 | L01;L02;L03;L04 | 1 | 4 | 83.16 |
| S_africanus | family2 | TB22_440_5_S6 | proband | Subsocial | DNBSEQ-G400 | V350235693 | L01;L02;L03;L04 | 1 | 4 | 83.15 |
| S_africanus | family2 | TB22_440_5_Fema | parent | Subsocial | DNBSEQ-G400 | V350235689 | L01;L02;L03;L04 | 1 | 4 | 83.1 |
| S_africanus | family2 | TB22_445_103_Ma | parent | Subsocial | DNBSEQ-G400 | V350235689 | L01;L02;L03;L04 | 1 | 4 | 82.9 |
| S_africanus | family3 | TB22_442_1_S1 | proband | Subsocial | DNBSEQ-G400 | V350236431 | L01;L02;L03;L04 | 1 | 4 | 83.16 |
| S_africanus | family3 | TB22_442_1_S2 | proband | Subsocial | DNBSEQ-G400 | V350236431 | L01;L02;L03;L04 | 1 | 4 | 80.15 |
| S_africanus | family3 | TB22_442_1_S3 | proband | Subsocial | DNBSEQ-G400 | V350231722 | L01;L02;L03;L04 | 1 | 4 | 68.82 |
| S_africanus | family3 | TB22_442_1_S4 | proband | Subsocial | DNBSEQ-G400 | V350231722 | L01;L02;L03;L04 | 1 | 4 | 83.12 |
| S_africanus | family3 | TB22_442_1_S5 | proband | Subsocial | DNBSEQ-G400 | V350231722 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family3 | TB22_442_1_S6 | proband | Subsocial | DNBSEQ-G400 | V350231722 | L01;L02;L03;L04 | 1 | 4 | 83.17 |
| S_africanus | family3 | TB22_442_1_Fema | parent | Subsocial | DNBSEQ-G400 | V350236431 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_africanus | family3 | TB22_442_7_Male | parent | Subsocial | DNBSEQ-G400 | V350236431 | L01;L02;L03;L04 | 1 | 4 | 83.16 |
| S_africanus | family4 | TB22_442_6_S1 | proband | Subsocial | DNBSEQ-G400 | V350223442 | L01;L02;L03;L04 | 1 | 4 | 83.11 |
| S_africanus | family4 | TB22_442_6_S2 | proband | Subsocial | DNBSEQ-G400 | V350223442;V350228912 | L01;L02;L03;L04 | 2 | 7 | 78.11 |
| S_africanus | family4 | TB22_442_6_S3 | proband | Subsocial | DNBSEQ-G400 | V350223445 | L01;L02;L03;L04 | 1 | 4 | 69.94 |
| S_africanus | family4 | TB22_442_6_S4 | proband | Subsocial | DNBSEQ-G400 | V350223445 | L01;L02;L03;L04 | 1 | 4 | 83.12 |
| S_africanus | family4 | TB22_442_6_S5 | proband | Subsocial | DNBSEQ-G400 | V350223445 | L01;L02;L03;L04 | 1 | 4 | 83.13 |
| S_africanus | family4 | TB22_442_6_S6 | proband | Subsocial | DNBSEQ-G400 | V350223445 | L01;L02;L03;L04 | 1 | 4 | 83.12 |
| S_africanus | family4 | TB22_442_4_Male | parent | Subsocial | DNBSEQ-G400 | V350223442 | L01;L02;L03;L04 | 1 | 4 | 83.15 |
| S_africanus | family4 | TB22_442_6_Fema | parent | Subsocial | DNBSEQ-G400 | V350223442 | L01;L02;L03;L04 | 1 | 4 | 83.14 |
| S_bicolor | family1 | S.bic_150.20_S1 | proband | Subsocial | DNBSEQ-G400 | V350166422 | L01 | 1 | 1 | 79.7 |
| S_bicolor | family1 | S.bic_150.20_S2 | proband | Subsocial | DNBSEQ-G400 | V350166422 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family1 | S.bic_150.20_S3 | proband | Subsocial | DNBSEQ-G400 | V350166422 | L03 | 1 | 1 | 81.07 |
| S_bicolor | family1 | S.bic_150.20_S4 | proband | Subsocial | DNBSEQ-G400 | V350166422 | L04 | 1 | 1 | 81.07 |
| S_bicolor | family1 | S.bic_150.20_S5 | proband | Subsocial | DNBSEQ-G400 | V350149154 | L01 | 1 | 1 | 81.08 |
| S_bicolor | family1 | S.bic_150.20_S6 | proband | Subsocial | DNBSEQ-G400 | V350149154 | L02 | 1 | 1 | 81.06 |
| S_bicolor | family1 | S.bic_150.20_female | parent | Subsocial | DNBSEQ-G400 | V350167939 | L04 | 1 | 1 | 81.06 |
| S_bicolor | family1 | S.bic_262.11_male | parent | Subsocial | DNBSEQ-G400 | V350149154 | L03 | 1 | 1 | 81.07 |
| S_bicolor | family2 | S.bic_150.6_S1 | proband | Subsocial | DNBSEQ-G400 | V350145631 | L04 | 1 | 1 | 81.09 |
| S_bicolor | family2 | S.bic_150.6_S2 | proband | Subsocial | DNBSEQ-G400 | V350145634 | L01 | 1 | 1 | 81.09 |
| S_bicolor | family2 | S.bic_150.6_S3 | proband | Subsocial | DNBSEQ-G400 | V350145634 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family2 | S.bic_150.6_S4 | proband | Subsocial | DNBSEQ-G400 | V350167939 | L01 | 1 | 1 | 81.07 |
| S_bicolor | family2 | S.bic_150.6_S5 | proband | Subsocial | DNBSEQ-G400 | V350167939 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family2 | S.bic_150.6_S6 | proband | Subsocial | DNBSEQ-G400 | V350167939 | L03 | 1 | 1 | 81.09 |
| S_bicolor | family2 | S.bic_150.6_female | parent | Subsocial | DNBSEQ-G400 | V350145631;V350168607 | L01;L02;L03;L04 | 2 | 5 | 72.72 |
| S_bicolor | family2 | S.bic_262.N5_male | parent | Subsocial | DNBSEQ-G400 | V350149154 | L04 | 1 | 1 | 81.07 |
| S_bicolor | family3 | 150.9_1 | proband | Subsocial | DNBSEQ-G400 | V350194601 | L03 | 1 | 1 | 81.07 |
| S_bicolor | family3 | 150.9_2 | proband | Subsocial | DNBSEQ-G400 | V350194601 | L04 | 1 | 1 | 81.08 |
| S_bicolor | family3 | 150.9_3 | proband | Subsocial | DNBSEQ-G400 | V350194542 | L01 | 1 | 1 | 81.08 |
| S_bicolor | family3 | 150.9_4 | proband | Subsocial | DNBSEQ-G400 | V350194661 | L04 | 1 | 1 | 81.08 |
| S_bicolor | family3 | 150.9_5 | proband | Subsocial | DNBSEQ-G400 | V350194704 | L01 | 1 | 1 | 81.09 |
| S_bicolor | family3 | 150.9_6 | proband | Subsocial | DNBSEQ-G400 | V350194704 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family3 | 150.9_female | parent | Subsocial | DNBSEQ-G400 | V350167206;V350168683 | L01;L04 | 2 | 2 | 85.42 |
| S_bicolor | family3 | 262.N3_male | parent | Subsocial | DNBSEQ-G400 | V350194632 | L04 | 1 | 1 | 81.09 |
| S_bicolor | family4 | 262.6_1 | proband | Subsocial | DNBSEQ-G400 | V350194550 | L01 | 1 | 1 | 81.08 |
| S_bicolor | family4 | 262.6_2 | proband | Subsocial | DNBSEQ-G400 | V350194550 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family4 | 262.6_3 | proband | Subsocial | DNBSEQ-G400 | V350194550 | L03 | 1 | 1 | 81.07 |
| S_bicolor | family4 | 262.6_4 | proband | Subsocial | DNBSEQ-G400 | V350194550 | L04 | 1 | 1 | 81.08 |
| S_bicolor | family4 | 262.6_5 | proband | Subsocial | DNBSEQ-G400 | V350194601 | L01 | 1 | 1 | 81.06 |
| S_bicolor | family4 | 262.6_6 | proband | Subsocial | DNBSEQ-G400 | V350194601 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family4 | 142.1_male | parent | Subsocial | DNBSEQ-G400 | V350194632 | L03 | 1 | 1 | 81.07 |
| S_bicolor | family4 | 262.6_female | parent | Subsocial | DNBSEQ-G400 | V350194597 | L01 | 1 | 1 | 81.08 |
| S_bicolor | family5 | 262.8_1 | proband | Subsocial | DNBSEQ-G400 | V350194597 | L03 | 1 | 1 | 81.09 |
| S_bicolor | family5 | 262.8_2 | proband | Subsocial | DNBSEQ-G400 | V350194597 | L04 | 1 | 1 | 81.09 |
| S_bicolor | family5 | 262.8_3 | proband | Subsocial | DNBSEQ-G400 | V350194555 | L01 | 1 | 1 | 81.09 |
| S_bicolor | family5 | 262.8_4 | proband | Subsocial | DNBSEQ-G400 | V350194555 | L02 | 1 | 1 | 81.07 |
| S_bicolor | family5 | 262.8_5 | proband | Subsocial | DNBSEQ-G400 | V350194555 | L03 | 1 | 1 | 81.06 |
| S_bicolor | family5 | 262.8_6 | proband | Subsocial | DNBSEQ-G400 | V350194555 | L04 | 1 | 1 | 81.08 |
| S_bicolor | family5 | 150.N8_male | parent | Subsocial | DNBSEQ-G400 | V350194632 | L02 | 1 | 1 | 81.08 |
| S_bicolor | family5 | 262.8_female | parent | Subsocial | DNBSEQ-G400 | V350194597 | L02 | 1 | 1 | 81.07 |
| S_lineatus | family1 | Slin-A8-S1 | proband | Subsocial | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 125.59 |
| S_lineatus | family1 | Slin-A8-S2 | proband | Subsocial | DNBSEQ-G400 | V300088301;V300088321;V300088397;V300088430 | L1;L2;L3;L4 | 4 | 7 | 96.22 |
| S_lineatus | family1 | Slin-A8-S3 | proband | Subsocial | DNBSEQ-G400 | V300088301;V300088321;V300088397;V300088430 | L1;L2;L3;L4 | 4 | 7 | 90.43 |
| S_lineatus | family1 | Slin-A8-S4 | proband | Subsocial | DNBSEQ-G400 | V300088301;V300088321;V300088397;V300088430 | L1;L2;L3;L4 | 4 | 7 | 89.83 |
| S_lineatus | family1 | Slin-A8-S5 | proband | Subsocial | DNBSEQ-G400 | V300088301;V300088321;V300088397;V300088430 | L1;L2;L3;L4 | 4 | 7 | 94.79 |
| S_lineatus | family1 | Slin-A8-S6 | proband | Subsocial | DNBSEQ-G400 | V300088301;V300088321;V300088397;V300088430 | L1;L2;L3;L4 | 4 | 7 | 96.57 |
| S_lineatus | family1 | Slin-A8-F | parent | Subsocial | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 113.69 |
| S_lineatus | family1 | Slin-E5-M | parent | Subsocial | DNBSEQ-G400 | V300088307;V300088321;V300088752;V300101732 | L1;L2;L3;L4 | 4 | 10 | 133.56 |
| S_lineatus | family2 | BOSK.16.S1 | proband | Subsocial | DNBSEQ-G400 | V350109501;V350109629 | L04 | 2 | 2 | 97.41 |
| S_lineatus | family2 | BOSK.16.S2 | proband | Subsocial | DNBSEQ-G400 | V350108243 | L01 | 1 | 1 | 90.87 |
| S_lineatus | family2 | BOSK.16.S3 | proband | Subsocial | DNBSEQ-G400 | V350108243 | L02 | 1 | 1 | 90.82 |
| S_lineatus | family2 | BOSK.16.S4 | proband | Subsocial | DNBSEQ-G400 | V350108243 | L03 | 1 | 1 | 90.79 |
| S_lineatus | family2 | BOSK.16.S5 | proband | Subsocial | DNBSEQ-G400 | V350108243 | L04 | 1 | 1 | 90.79 |
| S_lineatus | family2 | BOSK.16.S6 | proband | Subsocial | DNBSEQ-G400 | V350108233 | L01 | 1 | 1 | 91.08 |
| S_lineatus | family2 | BOSK.16.F | parent | Subsocial | DNBSEQ-G400 | V350108233 | L02 | 1 | 1 | 91.08 |
| S_lineatus | family2 | BOSK.16.M | parent | Subsocial | DNBSEQ-G400 | V350108233 | L03 | 1 | 1 | 91.13 |
| S_lineatus | family3 | BOSK.23.S1 | proband | Subsocial | DNBSEQ-G400 | V350109641 | L02 | 1 | 1 | 90.84 |
| S_lineatus | family3 | BOSK.23.S2 | proband | Subsocial | DNBSEQ-G400 | V350109641 | L03 | 1 | 1 | 90.84 |
| S_lineatus | family3 | BOSK.23.S3 | proband | Subsocial | DNBSEQ-G400 | V350109641 | L04 | 1 | 1 | 90.82 |
| S_lineatus | family3 | BOSK.23.S4 | proband | Subsocial | DNBSEQ-G400 | V350109501;V350109629 | L01 | 2 | 2 | 106.15 |
| S_lineatus | family3 | BOSK.23.S5 | proband | Subsocial | DNBSEQ-G400 | V350109501;V350109629 | L02 | 2 | 2 | 100.0 |
| S_lineatus | family3 | BOSK.23.S6 | proband | Subsocial | DNBSEQ-G400 | V350109501;V350109629 | L03 | 2 | 2 | 98.18 |
| S_lineatus | family3 | BOSK.23.F | parent | Subsocial | DNBSEQ-G400 | V350109703 | L01 | 1 | 1 | 91.13 |
| S_lineatus | family3 | BOSK.23.M | parent | Subsocial | DNBSEQ-G400 | V350109703 | L04 | 1 | 1 | 91.13 |
| S_lineatus | family4 | BOSK.10.S1 | proband | Subsocial | DNBSEQ-G400 | V350109668 | L04 | 1 | 1 | 90.82 |
| S_lineatus | family4 | BOSK.10.S2 | proband | Subsocial | DNBSEQ-G400 | V350109452 | L01 | 1 | 1 | 90.8 |
| S_lineatus | family4 | BOSK.10.S3 | proband | Subsocial | DNBSEQ-G400 | V350109452 | L02 | 1 | 1 | 90.85 |
| S_lineatus | family4 | BOSK.10.S4 | proband | Subsocial | DNBSEQ-G400 | V350109452 | L03 | 1 | 1 | 90.88 |
| S_lineatus | family4 | BOSK.10.S5 | proband | Subsocial | DNBSEQ-G400 | V350109452 | L04 | 1 | 1 | 90.8 |
| S_lineatus | family4 | BOSK.10.S6 | proband | Subsocial | DNBSEQ-G400 | V350109641 | L01 | 1 | 1 | 90.77 |
| S_lineatus | family4 | BOSK.10.F | parent | Subsocial | DNBSEQ-G400 | V350108207 | L03 | 1 | 1 | 91.14 |
| S_lineatus | family4 | BOSK.10.M | parent | Subsocial | DNBSEQ-G400 | V350108207 | L04 | 1 | 1 | 91.09 |
| S_lineatus | family5 | BOSK.13.S1 | proband | Subsocial | DNBSEQ-G400 | V350108233 | L04 | 1 | 1 | 91.11 |
| S_lineatus | family5 | BOSK.13.S2 | proband | Subsocial | DNBSEQ-G400 | V350103235;V350108206 | L01 | 2 | 2 | 91.14 |
| S_lineatus | family5 | BOSK.13.S3 | proband | Subsocial | DNBSEQ-G400 | V350103235;V350108206 | L02 | 2 | 2 | 91.05 |
| S_lineatus | family5 | BOSK.13.S4 | proband | Subsocial | DNBSEQ-G400 | V350103235;V350108206 | L03 | 2 | 2 | 91.06 |
| S_lineatus | family5 | BOSK.13.S5 | proband | Subsocial | DNBSEQ-G400 | V350103235;V350108206 | L04 | 2 | 2 | 91.04 |
| S_lineatus | family5 | BOSK.13.S6 | proband | Subsocial | DNBSEQ-G400 | V350108207 | L01 | 1 | 1 | 91.09 |
| S_lineatus | family5 | BOSK.13.F | parent | Subsocial | DNBSEQ-G400 | V350109668 | L03 | 1 | 1 | 90.85 |
| S_lineatus | family5 | BOSK.13.M | parent | Subsocial | DNBSEQ-G400 | V350108207 | L02 | 1 | 1 | 91.09 |
| S_tentoriicola | family1 | 91.S1 | proband | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 76.06 |
| S_tentoriicola | family1 | 92.S2 | proband | Subsocial | DNBSEQ-G400 | DP8400011525TR | L01 | 1 | 1 | 95.0 |
| S_tentoriicola | family1 | 93.S3 | proband | Subsocial | DNBSEQ-G400 | DP8400011525TR | L01 | 1 | 1 | 119.76 |
| S_tentoriicola | family1 | 94.S4 | proband | Subsocial | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 96.66 |
| S_tentoriicola | family1 | 95.S5 | proband | Subsocial | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 104.83 |
| S_tentoriicola | family1 | 96.S6 | proband | Subsocial | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 110.54 |
| S_tentoriicola | family1 | 81.F | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL;DP8400011737BL | L01 | 2 | 2 | 72.05 |
| S_tentoriicola | family1 | 86.M | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 78.32 |
| S_tentoriicola | family2 | 100.S4 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 99.97 |
| S_tentoriicola | family2 | 101.S5 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 102.31 |
| S_tentoriicola | family2 | 102.S6 | proband | Subsocial | DNBSEQ-G400 | DP8400011786TL | L01 | 1 | 1 | 109.7 |
| S_tentoriicola | family2 | 97.S1 | proband | Subsocial | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 103.82 |
| S_tentoriicola | family2 | 98.S2 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 87.85 |
| S_tentoriicola | family2 | 99.S3 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 86.68 |
| S_tentoriicola | family2 | 82.F | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL;DP8400011737BL | L01 | 2 | 2 | 73.14 |
| S_tentoriicola | family2 | 87.M | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 70.14 |
| S_tentoriicola | family3 | 103.S1 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 106.85 |
| S_tentoriicola | family3 | 104.S2 | proband | Subsocial | DNBSEQ-G400 | DP8400011678BR | L01 | 1 | 1 | 101.63 |
| S_tentoriicola | family3 | 105.S3 | proband | Subsocial | DNBSEQ-G400 | DP8400011786TL | L01 | 1 | 1 | 122.6 |
| S_tentoriicola | family3 | 106.S4 | proband | Subsocial | DNBSEQ-G400 | DP8400011454TL | L01 | 1 | 1 | 122.34 |
| S_tentoriicola | family3 | 107.S5 | proband | Subsocial | DNBSEQ-G400 | DP8400011516TR | L01 | 1 | 1 | 68.96 |
| S_tentoriicola | family3 | 108.S6 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 168.36 |
| S_tentoriicola | family3 | 83.F | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 67.63 |
| S_tentoriicola | family3 | 88.M | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 72.85 |
| S_tentoriicola | family4 | 109.S1 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 126.74 |
| S_tentoriicola | family4 | 110.S2 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 133.53 |
| S_tentoriicola | family4 | 111.S3 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 82.93 |
| S_tentoriicola | family4 | 112.S4 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 152.59 |
| S_tentoriicola | family4 | 113.S5 | proband | Subsocial | DNBSEQ-G400 | DP8400011786TL | L01 | 1 | 1 | 153.06 |
| S_tentoriicola | family4 | 114.S6 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 171.79 |
| S_tentoriicola | family4 | 84.F | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 65.67 |
| S_tentoriicola | family4 | 89.M | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 76.13 |
| S_tentoriicola | family5 | 115.S1 | proband | Subsocial | DNBSEQ-G400 | DP8400011679BR | L01 | 1 | 1 | 142.91 |
| S_tentoriicola | family5 | 116.S2 | proband | Subsocial | DNBSEQ-G400 | DP8400011376TR;FP100001127TL | L01 | 2 | 2 | 90.0 |
| S_tentoriicola | family5 | 117.S3 | proband | Subsocial | DNBSEQ-G400 | DP8400011516TR | L01 | 1 | 1 | 74.05 |
| S_tentoriicola | family5 | 118.S4 | proband | Subsocial | DNBSEQ-G400 | FP100001127TL | L01 | 1 | 1 | 68.57 |
| S_tentoriicola | family5 | 119.S5 | proband | Subsocial | DNBSEQ-G400 | FP100001127TL | L01 | 1 | 1 | 76.31 |
| S_tentoriicola | family5 | 120.S6 | proband | Subsocial | DNBSEQ-G400 | DP8400011516TR | L01 | 1 | 1 | 71.1 |
| S_tentoriicola | family5 | 85.F | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 66.06 |
| S_tentoriicola | family5 | 90.M | parent | Subsocial | DNBSEQ-G400 | DP8400011312BL | L01 | 1 | 1 | 75.21 |
