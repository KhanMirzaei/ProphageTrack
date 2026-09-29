# Caller setup and adapters

ProphageTrack does not install callers, download databases, upload sequences or connect to user accounts. Run callers in their own compatible environments. No default pair is selected until PhiSpy, PhageBoost and VIBRANT are benchmarked on relevant assemblies. geNomad is excluded.

## PhiSpy

Official documentation: https://github.com/linsalrob/PhiSpy

Requires an annotated GenBank genome, not bare FASTA. Provide matching FASTA and GenBank record IDs/sequences and real CDS annotations. Merely wrapping a FASTA in GenBank syntax is insufficient.

```bash
prophagetrack detect --caller phispy --sample-id S1 --assembly S1.fna \
  --genbank S1.gbk --out caller-results/S1_phispy --threads 2 \
  --caller-version YOUR_INSTALLED_VERSION
```

Or import `prophage_coordinates.tsv` directly:

```bash
prophagetrack import-calls --caller phispy --sample-id S1 \
  --input existing/prophage_coordinates.tsv --coordinates 1-inclusive --out S1_phispy.tsv
```

The adapter reads the first four fields (call ID, contig, start, end), preserving extra attachment information in the original caller output rather than claiming that ProphageTrack has evaluated it.

## PhageBoost

Official documentation: https://github.com/ku-cbd/PhageBoost

The current source and older packaged releases have different dependency requirements. Pin a tested release/commit and model hash in your study records. Use the region-level `phages_*.gff`, not the gene table `genecalls_*.gff3`.

Explicitly verify whether the installed version writes zero-based half-open or one-based inclusive endpoints, by comparing a known region and the exported sequence with the source assembly. Then choose the corresponding `--coordinates` value. This release does not guess from the filename.

```bash
prophagetrack import-calls --caller phageboost --sample-id S1 \
  --input existing/phages_S1.gff --coordinates 1-inclusive --out S1_phageboost.tsv
```

The example's `1-inclusive` is a declared convention, not a guarantee about every PhageBoost release. `detect` requires the same explicit option and supports `--executable /path/to/PhageBoost` for a separate environment.

## VIBRANT

Official documentation: https://github.com/AnantharamanLab/VIBRANT

VIBRANT is optional and its HMM searches have nontrivial resource costs. Its external environment and databases must already be configured. Use `integrated_prophage_coordinates.tsv`; free-virus predictions are not automatically integrated regions.

The TSV adapter defaults to column names `scaffold`, `start`, `stop`. If the installed release uses different names, specify them exactly:

```bash
prophagetrack import-calls --caller vibrant --sample-id S1 \
  --input existing/integrated_prophage_coordinates.tsv --coordinates 1-inclusive \
  --contig-column 'YOUR_CONTIG_HEADER' --start-column 'YOUR_START_HEADER' \
  --end-column 'YOUR_END_HEADER' --out S1_vibrant.tsv
```

`detect --caller vibrant` accepts these same column options, plus `--database`, `--caller-version` and `--database-version`. Its native output parsing has not been verified against a live run in this release.

## Failures versus no predictions

Some caller versions can exit successfully without creating expected output, including on configuration errors. ProphageTrack treats missing coordinate output as a failure requiring log inspection, not as an empty prediction set. A reviewed successful no-hit result can be represented by a header-only normalized TSV. Never convert an unexplained caller failure into zero prophages.

## Provenance and execution

Caller subprocesses receive argument lists, not shell command strings. They run with local staged input filenames. Their own internal shell behavior and dependencies remain external responsibilities. `caller_manifest.json` records the command, assembly hash, declared software/database versions and elapsed wall time. `caller.log` preserves output. Source files are not modified. Keep these manifests with analysis results; `run` does not infer missing caller/database versions from coordinate tables.

## Other tools

- MMseqs2: https://github.com/soedinglab/MMseqs2
- Pyrodigal: https://pyrodigal.readthedocs.io/
- Snakemake: https://snakemake.readthedocs.io/

Use upstream installation documentation and cite upstream methods used in the final study. Do not infer equivalent predictions across untested versions.
