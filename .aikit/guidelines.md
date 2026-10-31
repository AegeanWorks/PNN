# Project Guidelines

Rules for this project only, written by people or AI agents. They add to AIKit's principles and
rank below them: a guideline never overrides a principle (PRINCIPLE-000 R6, PRINCIPLE-001 R4).
Add one with `python aikit.py guideline <name> "<rule>" --source <user|ADR-NNN|person>`, or by
hand as one line below in the form `- **G-NNN** (source): rule`; IDs are never reused.

Persisted state carries a format version. A loader refuses a version it does not support,
without modifying or quarantining that file; older versions migrate when a migration exists.

## Rules
- **G-001** (ADR-001): Screening tiers are triage only; never report a compound as BBB-penetrant or efficacious without measured permeability and PNN-softening data.
- **G-002** (ADR-002): Docking consensus scores are relative to the docked set and never substitute for measured binding, PNN-softening or permeability data.
- **G-003** (ADR-003): Curated libraries keep source provenance; exact-SMILES deduplication is not chemical canonicalization and must be stated as such.
