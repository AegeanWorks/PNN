# Target Profiles & Receptor Invariants

## Purpose

`python -m pnn.target` provides programmatic access to structural biology invariants, binding-pocket coordinates, and AutoDock Vina grid parameter configurations for primary PNN-softening targets.

## Registered Receptor Targets

| Target ID | Gene | Primary PDB | Res (Å) | Pocket Type | Key Residues | Druggability |
|---|---|---|---|---|---|---|
| **`PTPsigma`** | `PTPRS` | 5K4P | 2.8 | `cspg_binding` | Lys201, Arg204, Lys207, Arg210, Tyr199, Phe203 | Moderate |
| **`PTPsigma_catalytic`** | `PTPRS` | 2PTS | 2.0 | `catalytic` | Cys267, Arg268, Asp326 | High |
| **`LAR`** | `PTPRF` | 2LAR | 1.9 | `catalytic` | Cys1522, Arg1528, Asp1490 | High |
| **`NgR1`** | `RTN4R` | 2OZ7 | 2.4 | `allosteric_lrr` | Arg68, Asn115, Arg140 | Moderate |
| **`MMP2`** | `MMP2` | 1CK7 | 1.5 | `catalytic` | His403, His407, His413, Glu404 | High |
| **`MMP9`** | `MMP9` | 1GKC | 2.3 | `catalytic_zinc` | His401, His405, His411, Glu402 | High |
| **`ADAMTS4`** | `ADAMTS4` | 2RJP | 1.5 | `catalytic_zinc` | His202, His206, His212, Glu203 | High |
| **`ChABC`** | `PL8` | 1L5D | 1.8 | `catalytic` | His200, Tyr247, Arg287 | Moderate |

## CLI Usage

```bash
# List all registered targets
python -m pnn.target list

# Display detailed receptor invariants and biological rationale
python -m pnn.target show PTPsigma

# Export AutoDock Vina grid configuration
python -m pnn.target config PTPsigma --receptor receptor.pdbqt --out docked.pdbqt --save vina.conf
```
