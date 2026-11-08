# E9 — parameter sensitivity of the family-level result

Question: does the headline family result (new genomes add about 1–3 points of novelty beyond a same-size baseline-vs-baseline null for
typical families, and about +10 points for a few enriched yeast families) depend on the clustering / region choices?

## Design

Five families, chosen to span the main cases:

| family | why | regions |
|---|---|---|
| `cazy_family_GH18` | CAZy, pilot family | 253,791 |
| `merops_family_S08` | MEROPS, pilot family | 264,977 |
| `pfam_PF00109` | secondary-metabolism backbone (ketosynthase) | 383,278 |
| `pfam_PF10450` | enriched Pfam (top enriched family with >= 5,000 regions: +0.10 at c70, interval [0.026, 0.183]) | 6,036 |
| `pfam_PF01568` | random-control Pfam (null difference about 0) | 3,964 |

Variants (all use the standard nested tiers c90 / c70 / c50 / c30; the tier-0 linclust is unchanged at 90 % coverage):

| variant | later-tier coverage | cov-mode | Pfam region |
|---|---|---|---|
| `env_cov80_m0` | 80 % | 0 (bidirectional) | envelope (default; reproduces the earlier wide run) |
| `env_cov60_m0` | 60 % | 0 | envelope |
| `env_cov90_m0` | 90 % | 0 | envelope |
| `env_cov80_m1` | 80 % | 1 (target-only) | envelope |
| `ali_cov80_m0` | 80 % | 0 | alignment span (`ali_from`..`ali_to`), Pfam families only |

The region definition applies only to Pfam (envelope vs alignment span). CAZy (dbCAN HMM `s_start`..`s_end`) and MEROPS (BLAST
`q_start`..`q_end`) have one span definition, so GH18 and S08 have no `ali` column (NA). The c30 tier is a default tier of
`family_cluster.py`; no extra run was needed to add it. Each variant re-ran F2 clustering and the F3 baseline-vs-baseline null
(`family_null.py`, 200 draws), then the table reports the new-minus-null difference in the novel fraction (genome-level, larger
matched size = 1,010 genomes per side) with the null-draw interval.

Code hooks (defaults unchanged): `family_cluster.py --cov-mode`, `family_extract.py --pfam-region {env,ali}`, launcher env vars
`COV`, `COVMODE`, `PFAMREGION`. Driver: `launch/e9_run.sh`; summary: `bin/e9_summary.py`; results: `results/e9/`
(`e9_summary.tsv`, per-variant `clusters_*`, `null_*`).

## Result: new-minus-null difference in the novel fraction (95 % interval)

| family | tier | env_cov80_m0 (default) | env_cov60_m0 | env_cov90_m0 | env_cov80_m1 | ali_cov80_m0 |
|---|---|---|---|---|---|---|
| GH18 | c70 | +0.010 [-0.001,+0.021] | +0.009 [-0.003,+0.019] | +0.009 [-0.003,+0.021] | +0.008 [-0.005,+0.019] | NA |
| GH18 | c50 | +0.002 [-0.002,+0.006] | +0.002 [-0.001,+0.005] | +0.004 [-0.000,+0.009] | +0.001 [-0.003,+0.005] | NA |
| GH18 | c30 | -0.000 [-0.001,+0.001] | +0.000 [-0.001,+0.001] | +0.002 [-0.000,+0.004] | +0.000 [-0.001,+0.001] | NA |
| S08 | c70 | +0.004 [-0.021,+0.025] | +0.000 [-0.023,+0.023] | +0.004 [-0.020,+0.025] | +0.001 [-0.024,+0.024] | NA |
| S08 | c50 | +0.000 [-0.014,+0.014] | +0.001 [-0.013,+0.014] | +0.001 [-0.013,+0.014] | +0.001 [-0.015,+0.013] | NA |
| S08 | c30 | -0.000 [-0.004,+0.004] | -0.001 [-0.004,+0.002] | -0.000 [-0.006,+0.005] | -0.001 [-0.006,+0.003] | NA |
| PF00109 | c70 | -0.000 [-0.006,+0.006] | +0.001 [-0.005,+0.006] | -0.001 [-0.007,+0.006] | +0.001 [-0.003,+0.006] | +0.001 [-0.005,+0.006] |
| PF00109 | c50 | -0.000 [-0.002,+0.002] | +0.000 [-0.001,+0.002] | -0.001 [-0.004,+0.002] | +0.000 [-0.001,+0.002] | +0.000 [-0.002,+0.002] |
| PF00109 | c30 | -0.000 [-0.002,+0.001] | -0.000 [-0.001,+0.000] | -0.001 [-0.003,+0.002] | +0.000 [-0.001,+0.001] | -0.000 [-0.001,+0.001] |
| PF01568 (control) | c70 | +0.001 [-0.036,+0.037] | +0.007 [-0.029,+0.040] | -0.002 [-0.041,+0.032] | +0.001 [-0.034,+0.040] | -0.005 [-0.040,+0.029] |
| PF01568 (control) | c50 | +0.001 [-0.019,+0.030] | +0.004 [+0.000,+0.038] | -0.000 [-0.028,+0.033] | +0.005 [-0.009,+0.039] | +0.008 [-0.002,+0.041] |
| PF01568 (control) | c30 | +0.001 [-0.009,+0.030] | +0.004 [+0.000,+0.038] | -0.002 [-0.026,+0.031] | +0.003 [-0.007,+0.037] | +0.003 [-0.006,+0.032] |
| PF10450 (enriched) | c70 | +0.100 [+0.026,+0.183] | +0.097 [+0.025,+0.180] | +0.101 [+0.032,+0.187] | +0.096 [+0.018,+0.187] | +0.095 [+0.025,+0.186] |
| PF10450 (enriched) | c50 | +0.071 [+0.010,+0.139] | +0.074 [+0.013,+0.133] | +0.074 [+0.008,+0.149] | +0.065 [-0.017,+0.130] | +0.070 [+0.010,+0.135] |
| PF10450 (enriched) | c30 | +0.018 [-0.018,+0.056] | +0.016 [-0.017,+0.050] | +0.018 [-0.029,+0.063] | +0.009 [-0.025,+0.042] | -0.007 [-0.043,+0.027] |

Counting rule used by `bin/e9_summary.py`: a *flip* = a variant whose difference has the opposite sign to the default (both beyond
0.003); a *size change* = a variant differing from the default by more than 0.01 and more than 50 % of the default difference.
Across 15 family x tier cells (each compared with 3 or 4 non-default variants) there is exactly one flip and one size change, both the
same cell: PF10450 at c30 with the alignment-span regions (-0.007 vs +0.018 for the default), where every interval spans 0.

## Interpretation

* The family-level result is **insensitive** to the coverage threshold (60 / 80 / 90 %), to bidirectional vs target-only coverage, and to
  envelope vs alignment-span Pfam regions. At c70 and c50 the variant estimates differ from the default by at most 0.004 for the
  large families, 0.005 for PF10450 and 0.007 for the small control, far inside the null-draw intervals.
* The typical-family pattern (differences of 0 to +1 point) holds in every variant for GH18, S08 and the PF00109 backbone family.
* The enriched family PF10450 keeps its excess (+0.095 to +0.101 at c70, intervals excluding 0 in all five variants; c50 +0.065
  to +0.074, the interval touches 0 only for target-only coverage). So the enrichment is not an artefact of one clustering setting.
* The 30 % tier adds no signal: all c30 differences are within +-0.02 with intervals spanning 0 (linclust/cluster is weak there,
  as recorded in the plan's open questions), so the tier set {90, 70, 50} is adequate for the headline.

## Caveats

* Five families only; PF10450 and PF01568 have 4–6 k regions, so their intervals are wide (control +-0.04 at c70).
* The default variant reruns the earlier wide-run settings and reproduces its numbers (GH18 +0.010, PF10450 +0.100 at c70).
* Only clustering and region definition were varied. Not varied: the mmseqs sensitivity (-s 6), the tier-0 coverage (90 %), the
  null design (K, draws) and the annotation thresholds applied upstream (Pfam GA, dbCAN, MEROPS cuts).
* Region definition could not be varied for CAZy and MEROPS (single span definition).
