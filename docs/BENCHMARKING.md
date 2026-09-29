# Benchmark before choosing default callers

The bundled synthetic demonstration does not measure prophage caller sensitivity or specificity. It contains fabricated prediction fixtures by design.

1. Select independently curated bacterial genomes representing the intended human-associated taxa, including negative regions and prophage remnants. Separate threshold-tuning genomes from held-out test lineages.
2. Run PhiSpy, PhageBoost and VIBRANT with fixed versions, databases and resource settings. Include annotation time/storage for PhiSpy and database setup/storage for all relevant callers.
3. Measure wall time and peak RAM separately for each caller and preprocessing stage, using `/usr/bin/time -v` on Linux or `/usr/bin/time -l` on macOS. Repeated runs and cache state should be reported. Detection manifests record wall time only; they do not claim peak-RAM measurement.
4. Normalize coordinates after verifying each tool's convention against actual sequence extraction.
5. Compare nucleotide precision/recall, region-level precision/recall, endpoint errors and tandem-prophage handling. The provided evaluator computes nucleotide overlap after interval union; region and boundary benchmarks require separate curation and analysis.
6. Repeat on controlled fragmentation and realistic incomplete assemblies. Evaluate missing predictions separately from unresolved loci.
7. Validate event calls using independently supported longitudinal cases. Include strain replacement, mixtures, repeated insertion sites, moved copies, Ns and missing flanks. Report unresolved fractions alongside error rates.
8. Choose the smallest caller set that delivers useful additional accuracy on held-out data. Do not equate majority voting with independent confirmation.

Example nucleotide benchmark:

```bash
prophagetrack evaluate --truth curated_regions.tsv --predictions predictions.tsv --out accuracy.json
```

Both inputs require `sample_id`, `contig`, `start`, `end` with 1-based inclusive coordinates. Truth must cover all evaluated samples, including explicit knowledge that unlisted samples/regions are negative. Incomplete truth can falsely penalize genuine predictions.

No real-data benchmark results, optimal caller ranking, resource guarantees or validated default thresholds are included in v0.1.0.
