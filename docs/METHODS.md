# Methods and interpretation

## Coordinates and caller evidence

All normalized calls, region tables and GFF outputs use 1-based inclusive coordinates. Python slices are converted internally. Caller imports require a declared convention; invalid coordinates, unknown contigs and duplicate call IDs within a sample/caller are rejected. Native PhageBoost coordinate semantics must be checked against the installed version; a .gff filename alone does not establish correctness.

Calls on the same sample/contig are linked if overlap covers at least 80% of each interval by default. A connected component is reconciled only if it is a complete graph and contains no repeated caller. Otherwise its calls remain separate with overlap flags. This avoids transitive chain merging and merging multiple same-caller regions. The representative interval is the observed call minimizing summed endpoint differences (ties broken by caller and ID). Neither union nor intersection invents a new boundary. All calls remain in `raw_calls.tsv`.

Different endpoints flag boundary disagreement. Any contig-edge contact flags possible truncation; edge contact does not prove incompleteness, and absence of edge contact does not prove completeness. Agreement counts are descriptive, not probabilities, and callers are not independent.

## Exact and related sequences

Unambiguous DNA is canonicalized against its reverse complement and SHA-256 labelled. This means 100% identity and full reciprocal length for extracted sequences. Ambiguous-base regions get occurrence-specific UNCERTAIN identifiers. Exact group labels do not prove complete prophage identity when boundaries are uncertain.

Exact sequence recovery searches all assemblies on both strands. It recovers caller misses, but alone cannot tell a prophage from a free viral contig, assembly contamination or another genomic context. Related nucleotide comparisons use MMseqs2 easy-search (nucleotide mode, both strands). Default thresholds are identity >=0.95 and reciprocal coverage >=0.85. If identity passes but only one sequence reaches coverage, the result is a fragment match. Matches covering neither sequence sufficiently are omitted from this pair table; shared genes can still connect them. The best local alignment is reported; disconnected rearranged blocks are not summed into whole-genome coverage. MMseqs2 heuristics may miss matches, especially short or divergent ones. There is no inferred phage species taxonomy or transitive related-group claim.

## Genes and protein families

Pyrodigal in metagenomic mode predicts genes consistently in extracted regions; partial-gene flags are retained. Alternatively, supplied proteins require a one-to-one mapping to region IDs and relative 1-based inclusive coordinates. Supplied annotation content and translation correctness remain the provider's responsibility; coordinate bounds and sequence IDs are checked.

Protein families are greedy representative groups. Members must match the representative at default identity >=0.50 and reciprocal coverage >=0.80 using MMseqs2, or have identical protein sequences. Groups are not orthology assertions and not necessarily pairwise cliques. Exact mode groups only identical proteins. Thresholds are exploratory and should be tested on the target taxa. Three consecutive family IDs shared in either order are reported as conserved windows; overlapping windows are not merged into reconstructed recombination events. Boundary contamination, common mobile-element proteins and annotation quality can influence sharing.

## Tracking

A cohort is a subject/strain pair whose samples have dates and nonempty lineage_evidence. The software trusts this supplied assignment; it does not determine strain relatedness. Evidence should come from a suitable independent core-genome/SNP analysis, not prophage sharing itself or species-level ANI alone. Within-host strain mixtures and replacement can invalidate naive longitudinal interpretation.

Only unflagged reference regions with sufficient flanking sequence become anchors. Default flank length is 200 bp on each side. Both flanks must contain only A/C/G/T and each must match exactly once across both strands of the target assembly. The two hits must be collinear on the same contig. Missing, repeated, differently oriented or split flanks produce unresolved status. An interval with ambiguous bases is unresolved.

An intervening sequence identical to the reference is occupied. A short resolved gap (default <=50 bp) is candidate_empty, allowing boundary/attachment-site remnants. This does not assert molecular excision and may miss substantial remnants; set --max-empty-gap 0 for a strict abutting-flank criterion. A longer differing interval is different_occupancy, not automatic absence. Reference intervals at or below the empty-gap limit are ineligible.

Consecutive eligible samples are classified as follows:

- occupied -> occupied: persistence at the corresponding site;
- candidate_empty -> occupied: candidate_gain;
- occupied -> candidate_empty: candidate_loss;
- occupied <-> different_occupancy: candidate_modification_or_replacement;
- candidate_empty -> candidate_empty: no_tracked_prophage_at_site;
- other transitions: unresolved.

Same-date replicates are not ordered into biological events. Loci sharing an exact flank pair, reference sequence and cohort are deduplicated; divergent reference boundaries can yield distinct loci, so event counts must not be interpreted as automatically independent events. An unknown intermediate sample is not skipped to infer an event across it. A local empty site does not rule out the prophage elsewhere in the same assembly; whole-assembly recovery is a separate table.

No read mapping, strain inference, induction assessment, host-range prediction, transmission direction or causal analysis is implemented in v0.1.0. These require separate evidence and validation. Results describe candidate locus changes and observed sharing.

## Scaling

Exact recovery scans each distinct sequence against each assembly. Related comparisons are all-versus-all across extracted regions. These implementations are suitable for pilot studies; no population-scale RAM/runtime guarantee is established. Use a representative pilot before large cohorts and report limits separately for annotation, each caller and downstream comparisons. No hard-coded all-callers default is claimed to be optimal.
