# Validation status — 0.1.0 research preview

Validated locally on 2026-09-29. This is software validation, not a biological validation claim.

## Completed

- 39 automated tests passed using Python 3.11.4 and Biopython 1.87.
- Tests cover coordinate bounds, unknown/duplicate inputs, reverse-complement identity, ambiguous DNA, caller boundary disagreements, transitive overlap chains, repeated/split/missing flanks, ambiguous gaps, changed occupancy, exact recovery in caller-negative samples, gene sharing, HTML escaping, empty predictions and refused overwrites.
- Synthetic longitudinal analysis returned the expected persistence, candidate loss and candidate gain, with no automatic event calls when lineage evidence was missing.
- Synthetic imported proteins produced five exact families across three regions, including three families shared with the different prophage sequence.
- MMseqs2 18-8cc5c nucleotide and protein searches executed successfully on the demonstration; protein-family results matched the intended shared genes. The different prophage did not pass the default whole-sequence coverage threshold, as expected from its shared three-of-four-gene design.
- Pyrodigal 3.7.1 in an existing local Python environment executed successfully on synthetic regions and produced 12 predicted genes and five exact protein families.
- Three external-caller runner tests used mock executables to verify argument handling, provenance and rejection of exit-zero runs with no output. These are not live PhiSpy, PhageBoost or VIBRANT tests.
- Source distribution and wheel built locally. The wheel installed without network downloads and its CLI completed the synthetic workflow outside the source directory.
- Snakemake 9.23.1 dry-run passed and built the expected two-job DAG. The workflow itself was not executed end-to-end with external callers.

## Not yet validated / not implemented

- No live execution of PhiSpy, PhageBoost or VIBRANT was performed. Installation/version compatibility, real native output parsing and no-hit behavior need version-specific integration tests.
- No real bacterial cohort, independently curated prophage truth or resource comparison was evaluated. No best caller pair is claimed.
- No induction, infection, host range, horizontal-transfer direction, read-based evidence or automatic bacterial strain assignment is inferred.
- Exact-flank tracking will miss some genuine events in divergent or fragmented assemblies. Threshold sensitivity and unresolved rates need real-data assessment.
- This release provides one analysis-level Snakemake rule, not full internal checkpointing.

