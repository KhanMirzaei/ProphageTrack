# ProphageTrack

**Identify shared prophage sequences and genes, and assess prophage persistence, candidate gain and candidate loss within human-associated bacterial strains.**

Version 0.1.0 is a **research preview**. The comparative and tracking code is functional and tested with synthetic cases. PhiSpy, PhageBoost and VIBRANT adapters are provided, but their live execution and biological accuracy have not yet been validated for this release. No default caller pair has been benchmark-selected. geNomad is not used.

## What it does

- Imports prophage coordinates from PhiSpy, PhageBoost, VIBRANT or a normalized TSV.
- Can launch each locally installed caller explicitly, preserving logs and provenance.
- Retains each caller's predictions and reconciles only unambiguous overlapping groups.
- Groups identical extracted DNA sequences, including reverse complements. Regions containing ambiguous bases remain separate.
- Uses optional MMseqs2 nucleotide alignments to report related sequences and fragment matches.
- Searches every assembly for full exact copies of discovered sequences, including caller-negative samples.
- Predicts genes with Pyrodigal, or accepts supplied proteins; identifies exact proteins or homologous protein families with MMseqs2 and shared three-gene windows.
- Projects unique bacterial flanks across a user-supported longitudinal strain lineage to assess persistence, candidate gain/loss and changed occupancy.
- Produces TSV matrices, sequences, GFF3, network edges, evidence tables and an offline HTML report.

## Scope and limitations

An identical extracted region is not automatically a complete prophage. "Not detected" is not verified absence. Sequence sharing is not proof of horizontal transfer, infection, induction or direction of transmission. The tool does not infer strain identity: a strain name plus meaningful independently established `lineage_evidence` is required for tracking. It does not currently analyze reads, excision junctions, viromes, or host infectivity.

Exact flank matching is deliberately conservative: mutations, repeats and assembly gaps can make sites unresolved. Boundary-disagreeing, overlapping, ambiguous and edge-truncated regions are excluded as tracking anchors. Related prophages are reported as pairwise alignments, not forced into transitive identity clusters. Search thresholds are exploratory, not species definitions.

## Install

Python 3.10 or later is required. From a downloaded copy of this repository:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install '.[genes]'
```

Install MMseqs2 separately for nucleotide and protein comparisons. External callers each require their own compatible installation and, where applicable, databases. Nothing is downloaded automatically by ProphageTrack. See [caller setup](docs/CALLERS.md).

For an analysis using supplied proteins and exact comparisons only, `pip install .` is sufficient. Run `prophagetrack --help` for commands. Linux is the intended platform for the complete external-caller workflow; the analysis core has also been tested on macOS. Windows users should use a Linux environment such as WSL.

## Try the synthetic demonstration

The repository includes small, completely synthetic assemblies, predicted regions and proteins in `examples/synthetic/`. They test analysis logic, not caller performance.

```bash
prophagetrack run \
  --samples examples/synthetic/samples.tsv \
  --calls examples/synthetic/calls.tsv \
  --proteins examples/synthetic/proteins.faa \
  --gene-map examples/synthetic/gene_map.tsv \
  --related none --gene-mode exact \
  --label 'Entirely synthetic demonstration' \
  --out results/demo
```

Open `results/demo/report.html`. Expected outcomes are persistence from A0 to A1, candidate loss from A1 to A2, and candidate gain from A2 to A3. A3 has no prediction in the supplied caller table: the assembly search and site comparison recover it. B0 contains a different sequence sharing three of four synthetic genes. The fictional strain assignments are fixture metadata, not biological inference.

With MMseqs2 installed, remove `--related none --gene-mode exact` to run alignments. To test real gene prediction, also remove `--proteins` and `--gene-map` after installing the `genes` extra.

## Your own data

Create a tab-separated sample sheet:

```text
sample_id	assembly	subject_id	strain_id	date	lineage_evidence
person1_day0	assemblies/day0.fna	person1	strain_A	2025-01-01	independent_core_genome_analysis
person1_day30	assemblies/day30.fna	person1	strain_A	2025-01-31	independent_core_genome_analysis
```

The displayed `\t` separators above represent real tab characters. Assembly paths are relative to the sample sheet. Each row is one assembly/genome or MAG; `sample_id` must uniquely identify that row. Cross-sectional inputs need only `sample_id` and `assembly`. Subject, strain, ISO date and lineage evidence are all required to include a sample in tracking. Strain IDs are scoped within subjects; matching names across subjects do not establish the same lineage. Unassigned contigs can be compared but should not be given invented strain labels.

Normalize caller outputs with `import-calls`, or create this schema using **1-based inclusive** coordinates:

```text
sample_id	contig	start	end	caller	call_id
person1_day0	contig_1	12001	48500	phispy	pp1
```

Then run:

```bash
prophagetrack run --samples samples.tsv --calls phispy_calls.tsv phageboost_calls.tsv \
  --out results/study --threads 4
```

Use `--skip-genes` explicitly if gene analysis is not wanted. Output directories must be new or empty to prevent mixing runs. An interrupted run should be moved aside before restarting. `COMPLETE` is written only after successful completion. The optional Snakemake wrapper resumes at the completed analysis-rule level; this release does not checkpoint individual internal analysis stages.

## Outputs

| File | Contents |
|---|---|
| `raw_calls.tsv`, `regions.tsv` | Original calls, chosen intervals, caller support and uncertainty flags |
| `prophages.fna`, `prophages.gff3` | Extracted sequences and assembly coordinates |
| `genome_prophage_matrix.tsv` | Called, recovered or not-detected sequences; uncertainty is explicit |
| `exact_recovery.tsv`, `related_pairs.tsv` | Full exact assembly matches and pairwise alignment evidence |
| `group_prevalence.tsv` | Counts of assemblies, named within-subject strains and subjects; not population prevalence estimates |
| `site_observations.tsv`, `events.tsv` | Projected loci and candidate longitudinal transitions |
| `genes.tsv`, `proteins.faa`, `gene_families.tsv` | Gene coordinates and families |
| `prophage_gene_matrix.tsv`, `shared_gene_blocks.tsv` | Gene copy counts and conserved family windows |
| `network_edges.tsv` | Sample–region–sequence group–gene family links for Cytoscape or other software |
| `manifest.json`, `report.html` | Parameters, input hashes, versions, limitations and readable results |

Matrices and comparisons refer to extracted regions and include uncertain predictions with flags. Evidence tables must be consulted before biological interpretation. Reports and databases are generated outputs and should not be committed to the repository.

## Validation and reproducibility

Run `python -m unittest discover -s tests -v`. See [methods](docs/METHODS.md), [validation status](docs/VALIDATION.md), and [benchmark plan](docs/BENCHMARKING.md). Exact dependency versions used for an analysis are recorded in `manifest.json`; external-caller version and database details must also be supplied and retained in the detection manifests.

For the optional workflow wrapper, install `.[genes,workflow]`, edit `workflow/config.yaml`, then run:

```bash
snakemake --snakefile workflow/Snakefile --cores 4
```

## Repository contents

This repository contains source code, tests, installation metadata, workflow configuration, scientific documentation and a small synthetic example. The example is retained so users can verify installation without downloading a biological dataset.

Generated HTML reports, analysis results, databases, local environments and distribution archives belong outside the repository. The report renderer in `prophagetrack/report.py` is source code and is required to generate reports.

## Attribution

ProphageTrack integrates independently developed tools; it does not replace their methods or citations. Please cite every caller and comparison tool used, in addition to ProphageTrack. Links are provided in [CALLERS.md](docs/CALLERS.md). ProphageTrack's own code is MIT licensed; external tools and databases retain their own licenses.
