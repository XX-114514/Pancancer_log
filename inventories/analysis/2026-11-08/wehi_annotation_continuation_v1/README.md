# New TME: verified R3 continuation

[Method](../../../../methods/TME-WEHI-R3-consensus-continuation-20261108.md) · [Run](../../../../runs/20261108_WEHI_R3_continuation_v1.md)

This is an evidence-based working annotation, not a final biological release. All matrices and per-cell products remain on the research server.

| Pool | Seed-file cells | Stable at resolution 1 | Main subtype candidates |
|---|---:|---:|---:|
| T_NK | 1,415,387 | 1,112,503 | 811,800 |
| B_Plasma | 297,478 | 238,504 | 95,472 |
| Stromal | 243,620 | 183,771 | 101,013 |
| Myeloid | 656,808 | 452,373 | 328,756 |
| Endothelial | 137,014 | 99,893 | 50,895 |

Main candidates total 1,387,936. These are neither all TME cells nor the CoVarNet eligibility denominator. The B seed object lacks 2,384 cells reassigned from Myeloid; they retain explicit pending subtype status. All five full seed-label mappings, exchange fractions and stability values passed a separate pandas-crosstab replay.

| New transfer | Query cells | Unresolved | Holdout macro accuracy |
|---|---:|---:|---:|
| Microglia | 12,361 | 997 | 0.945 |
| Endothelial | 32,088 | 10,750 | 0.845 |
| Myeloid | 84,617 | 18,373 | 0.871 |

Holdout accuracy is an upper-bound diagnostic, not real-query accuracy. Transferred labels remain sensitivity-only. The literal consensus merge rule produced mixed identity/state/sample-specific components; these are kept in diagnostics but excluded from main subtype candidates.

All 41 R3 clusters have working actions, with 388 Myeloid and 1,016 EC cells set aside. The current sinusoidal-like EC signal is cluster12; the old task summary referred to cluster6. Rare T/NK splits remain observed candidates because resolution-1 stability does not validate their finer resolution. Four tissue pools lack ten-seed consensus. October comparison and downstream scope/threshold acceptance remain open.

The full 3,347,659-cell working sidecar passed set-preservation, exact main-candidate membership and status readback at 06:22:10 UTC. Six pools total 535,749 transferred and 122,572 unresolved low-depth cells; both categories are excluded from main proportions. See [acceptance](working_acceptance.json).
