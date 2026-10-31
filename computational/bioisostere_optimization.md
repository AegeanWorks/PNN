# Bioisosteric Replacement & Metabolic Soft-Spot Optimization

## Overview

The `pnn.bioisostere` module applies a fixed table of medicinal-chemistry transformation rules to
candidate SMILES strings, recomputes continuous Pfizer CNS MPO scores from rule-estimated property
deltas, and flags optimization potential.

Matching is **literal substring replacement**, not chemical canonicalization. A rule is skipped
when the matched motif would be cut out of a ring context or the result would leave unbalanced ring
labels. Molecular weight, logP, TPSA, HBD, HBA and rotatable-bond values are rule-estimated deltas,
not recomputed from the analog structure; every analog requires cheminformatics validation and
synthesis before any potency, stability or exposure claim.

## Codified Transformation Rules

| Rule | Target motif | Replacement | Rationale |
|---|---|---|---|
| `BIO-01` | `c1ccccc1` | `c1ccc(F)cc1` | Para-fluorination blocks CYP aromatic hydroxylation. |
| `BIO-02` | `c1ccccc1` | `c1ccc(C(F)(F)F)cc1` | CF3 confers metabolic inertness and higher lipophilicity. |
| `BIO-03` | `NC(=O)C` | `NC(=O)C1CC1` | Cyclopropyl constrains the methyl group and rigidifies the amide. |
| `BIO-04` | `N1CCNCC1` | `N1CCOCC1` | Morpholine oxygen lowers the distal amine pKa by ~1.5 units. |
| `BIO-05` | `CN1CCNCC1` | `C([2H])([2H])N1CCNCC1` | Benzylic dideuteration slows CYP N-dealkylation (kinetic isotope effect). |
| `BIO-06` | `C(=O)OC` | `c1nc(C)no1` | 1,2,4-Oxadiazole replaces an esterase-labile ester. |

Ring labels in a motif are relabelled to match the parent's numbering, and replacement ring labels
are renumbered to a free digit. This keeps generated SMILES ring-balanced; it does not make the
substitution chemically canonical.

## Property Deltas

Each rule carries an explicit delta row (`delta_mw`, `delta_logp`, `delta_tpsa`, `delta_hbd`,
`delta_hba`, `delta_rotb`). Examples:

- `BIO-03` adds 40.0 Da, +0.60 logP, removes one rotatable bond and does not change HBA.
- `BIO-04` adds 1.0 Da, −0.35 logP, −2.8 Å² TPSA, −1 HBD and no HBA change.
- `BIO-05` adds 2.0 Da with no other change.

The CNS MPO score is recomputed from these estimates, so it ranks hypotheses within the rule table
rather than predicting measured properties.

## Skipped Transformations

A rule does not fire when the motif is immediately followed by a ring label: replacing that atom
would splice the label onto an unrelated replacement atom and corrupt the ring. For example, in
`...NC(=O)C2=CC=C...` the amide carbon is part of the aromatic ring, so `BIO-03` is skipped instead
of producing an invalid analog.

## Command-Line Usage

```bash
python -m pnn.bioisostere \
  --compounds data/compounds/vendor_library.csv \
  --parent PNN-101 \
  --out output/bioisosteric_analogs.csv \
  --report output/bioisostere_report.md \
  --json output/bioisostere_summary.json
```

The analog CSV records `rule_id`, both SMILES strings and the estimated deltas. The Markdown report
repeats the method limits so generated analogs are never presented as validated molecules.
