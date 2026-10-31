# Compound Library Curation

## Purpose

`python -m pnn.curate` merges vendor, repurposing and in-house compound CSVs into
one traceable screening library. It screens every row with the property rules,
removes duplicate identities, applies the tier floor, and caps how many compounds
each region of property space may contribute.

## Inputs

Each input is a UTF-8 CSV using the screening contract
([screening_pipeline.md](screening_pipeline.md)): the required columns are
`compound_id`, `mw`, `logp`, `tpsa`, `hbd`, `hba` and `rotatable_bonds`; `smiles`,
`name` and the optional ADMET fields are carried when present. Unlike the screening
command, curation accepts duplicate rows inside one file and resolves them.

## Deduplication

1. **By `compound_id`**: rows with the same ID are collapsed to the better screening
   result (tier first, then score). If two non-empty SMILES disagree for one ID,
   curation refuses the merge because the identity is ambiguous.
2. **By exact SMILES**: after the ID pass, compounds with the same non-empty SMILES
   string are collapsed to the better result. Blank SMILES are never used to merge
   different IDs.

Exact string matching is **not** chemical canonicalization: tautomers, salts,
stereochemistry variants and differently written SMILES remain separate. Use a
chemistry toolkit to canonicalize before curation when that matters.

## Tier Filter

`--min-tier` (default `B`) keeps that screening tier and better. `Reject` rows are
never emitted at any setting.

## Property-Space Diversity

Six core properties are normalized against the screening windows:

| Property | Range used for normalization |
|---|---|
| `mw` | 200-450 Da |
| `logp` | 1.5-3.5 |
| `tpsa` | 0-70 A2 |
| `hbd` | 0-3 |
| `hba` | 0-7 |
| `rotatable_bonds` | 0-7 |

Values outside a range are clamped to its edge. Each normalized property is split
into `--bin-width` cells (default `0.1`, i.e. ten cells per axis), giving a
six-dimensional bin label. At most `--per-bin` compounds (default `25`) survive
from each bin, selected by tier, then score, then ID. `--max-compounds` applies a
final global cap. This is property-sphere coverage, not scaffold or chemotype
diversity; it keeps the set from clustering around one property profile.

## Usage

```bash
python -m pnn.curate vendor_a.csv vendor_b.csv \
    --out curated_library.csv --json curation_summary.json \
    --min-tier B --per-bin 25 --bin-width 0.1 --max-compounds 200000
```

- `--out` defaults to `curated_library.csv`; parent directories are created.
- `--json` writes input counts, duplicate counts, tier exclusions, bin count,
  selected rows and whether the global cap truncated the set.
- Exit codes: `0` success, `1` invalid input/conflict/output, `2` misuse.

## Output

`compound_id`, `name`, `smiles`, `source`, `tier`, `score`, `mw`, `logp`, `tpsa`,
`hbd`, `hba`, `rotatable_bonds`, `pka`, `charge`, `property_bin`,
`failed_criteria`, `unknown_criteria`.

`source` is the input path each kept compound came from; `property_bin` is the
`-`-joined six-cell label. The JSON summary is the audit trail for every removal
count.

## Limits

- Curation prioritizes by the Phase 2 property heuristics. It does not establish
  activity, permeability or safety; curated compounds still need docking and the
  Phase 3 assays.
- Exact-SMILES deduplication can miss salt forms and tautomers; canonicalize
  upstream when duplicate chemistry matters.
- Per-bin caps are an approximation of property-space spacing, not a maximum
  distance guarantee.
- The kept record's provenance is one source, even when the same compound appeared
  in several libraries; the JSON duplicate counts record the collapse.
- A fixed input set and fixed options are required to reproduce a library.
