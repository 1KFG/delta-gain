#!/usr/bin/env nextflow
/*
 * protein_novelty.nf — DeltaGain Stage 3.
 *
 * Three pieces, in dependency order:
 *   1. FETCH_UNIPROT_FUNGI (REST /stream, reviewed-only) — small connectivity/
 *      content smoke test. Live-tested 2026-09-06.
 *   2. FETCH_UNIPROT_FUNGI_FTP -> CONVERT_DAT_TO_FASTA -> COMBINE_UNIPROT_FASTA
 *      -> MMSEQS_LINCLUST -> DIAMOND_MAKEDB — the real, comprehensive
 *      UniProtKB-fungi (Swiss-Prot + TrEMBL, ~17.5M sequences) reference set,
 *      bulk-downloaded from FTP (not REST -- see FETCH_UNIPROT_FUNGI_FTP's
 *      header comment for why /stream can't serve this), converted to FASTA
 *      + InterPro/Pfam/GO cross-refs, combined, dereplicated via MMseqs2
 *      linclust, then built into a DIAMOND db for step 3.
 *   3. BUILD_PEP_MANIFEST -> COMBINE_BFD_PROTEINS -> MMSEQS_LINCLUST (again,
 *      BFD's own proteins this time) -> SPLIT_FASTA -> [DIAMOND_NOVELTY_SCREEN
 *      vs. nr_cluster_seq.dmnd, DIAMOND_SEARCH_UNIPROT vs. step 2's output]
 *      -> BIN_NOVELTY_HITS — the BFD-side novelty screen (HANDOFF.md items
 *      1/3/4): extract BFD's annotated proteins (input/pep/, mapped to ASMID
 *      via BUILD_PEP_MANIFEST -- see its header for why the pep filenames
 *      aren't already ASMID-keyed), dereplicate them the same way the
 *      UniProt-fungi set was, then screen the non-redundant representatives
 *      against both reference sets and bin the hits (E-value + coverage,
 *      per the design doc's framing fix).
 *
 * Step 2 update (2026-09-08): it WAS run for real (before the DeltaGain/
 * DeltaGain_Fungi repo split) -- 15,089,757 representative sequences, at
 * `DeltaGain/results/protein_novelty/uniprot_fungi_nr.fasta` (the OLD
 * pre-split location, not DeltaGain_Fungi's `results/`). Re-deriving it from
 * scratch (13 GB FTP download + ~17.5M-sequence linclust, 48h) would be pure
 * waste when that file already exists, so `params.uniprot_fungi_nr_fasta`
 * (set in profile_protein_novelty.config) short-circuits step 2 straight to
 * DIAMOND_MAKEDB when non-empty. Set it to '' explicitly to force a real
 * re-derivation (e.g. picking up a newer UniProt release).
 *
 * Step 3 is sized for a pilot run first (params.n_test, see
 * profile_protein_novelty.config's safety default) against the several
 * thousand genomes currently annotated in Fungi_BFD_runs/input/pep/, not the
 * full ~22k set (most of which isn't annotated yet -- BUILD_PEP_MANIFEST
 * reports exactly how many are, each run).
 *
 * DIAMOND_CLUSTER + COMPARE_CLUSTERINGS (added 2026-09-08): a second,
 * parallel dereplication of the SAME combined BFD protein set, using
 * DIAMOND's cascaded clustering (the algorithm behind DIAMOND DeepClust,
 * Buchfink et al. 2026 Nat Methods) instead of MMSEQS_LINCLUST_BFD's
 * MMseqs2 linclust -- run side by side, NOT wired into the search/binning
 * path, purely to answer "does this cluster BFD's proteins meaningfully
 * differently?" via COMPARE_CLUSTERINGS' Adjusted Rand Index. See
 * modules/dedup/DIAMOND_CLUSTER's header for a real DIAMOND 2.2.6 bug
 * found and worked around while building this (`--reps`+`--out` crashes).
 */

nextflow.enable.dsl = 2

include { FETCH_UNIPROT_FUNGI }     from './modules/refdb/FETCH_UNIPROT_FUNGI/main.nf'
include { FETCH_UNIPROT_FUNGI_FTP } from './modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf'
include { CONVERT_DAT_TO_FASTA }    from './modules/refdb/CONVERT_DAT_TO_FASTA/main.nf'
include { COMBINE_UNIPROT_FASTA }   from './modules/refdb/COMBINE_UNIPROT_FASTA/main.nf'
include { MMSEQS_LINCLUST as MMSEQS_LINCLUST_UNIPROT } from './modules/dedup/MMSEQS_LINCLUST/main.nf'
include { MMSEQS_LINCLUST as MMSEQS_LINCLUST_BFD }     from './modules/dedup/MMSEQS_LINCLUST/main.nf'
include { DIAMOND_MAKEDB }          from './modules/search/DIAMOND_MAKEDB/main.nf'
include { BUILD_PEP_MANIFEST }      from './modules/pep/BUILD_PEP_MANIFEST/main.nf'
include { COMBINE_BFD_PROTEINS }    from './modules/pep/COMBINE_BFD_PROTEINS/main.nf'
include { SPLIT_FASTA }             from './modules/split/SPLIT_FASTA/main.nf'
include { DIAMOND_NOVELTY_SCREEN }  from './modules/search/DIAMOND_NOVELTY_SCREEN/main.nf'
include { DIAMOND_SEARCH_UNIPROT }  from './modules/search/DIAMOND_SEARCH_UNIPROT/main.nf'
include { BIN_NOVELTY_HITS }        from './modules/binning/BIN_NOVELTY_HITS/main.nf'
include { DIAMOND_CLUSTER; EXTRACT_CLUSTER_REPS } from './modules/dedup/DIAMOND_CLUSTER/main.nf'
include { COMPARE_CLUSTERINGS }     from './modules/binning/COMPARE_CLUSTERINGS/main.nf'
include { BUILD_CLUSTER_INCIDENCE_MATRIX } from './modules/binning/BUILD_CLUSTER_INCIDENCE_MATRIX/main.nf'
include { ACCUMULATION_CURVE }             from './modules/binning/ACCUMULATION_CURVE/main.nf'

workflow {
    // Small connectivity/content smoke test (REST) -- cheap, safe at any scale.
    FETCH_UNIPROT_FUNGI(Channel.of('reviewed'))

    // Step 2: the real comprehensive UniProt-fungi reference set. Reuse the
    // already-built real result by default (see header comment) rather than
    // re-deriving it -- only runs the FTP/convert/combine/linclust chain
    // when params.uniprot_fungi_nr_fasta is explicitly cleared.
    if (params.uniprot_fungi_nr_fasta) {
        uniprot_nr_ch = Channel.fromPath(params.uniprot_fungi_nr_fasta, checkIfExists: true)
    } else {
        divisions_ch = Channel.of(
            tuple('sprot',  'uniprot_sprot_fungi.dat.gz'),
            tuple('trembl', 'uniprot_trembl_fungi.dat.gz'),
        )

        FETCH_UNIPROT_FUNGI_FTP(divisions_ch)
        CONVERT_DAT_TO_FASTA(FETCH_UNIPROT_FUNGI_FTP.out.dat_gz)

        COMBINE_UNIPROT_FASTA(
            CONVERT_DAT_TO_FASTA.out.map { division, fasta, xrefs -> fasta }.collect(),
            CONVERT_DAT_TO_FASTA.out.map { division, fasta, xrefs -> xrefs }.collect(),
        )

        MMSEQS_LINCLUST_UNIPROT(COMBINE_UNIPROT_FASTA.out.fasta.map { fasta -> tuple('uniprot_fungi', fasta) })
        uniprot_nr_ch = MMSEQS_LINCLUST_UNIPROT.out.representatives.map { prefix, fasta -> fasta }
    }

    DIAMOND_MAKEDB(uniprot_nr_ch)

    // Step 3: BFD's own annotated proteins -- extract, dereplicate, screen.
    samples_csv = file(params.samples, checkIfExists: true)

    BUILD_PEP_MANIFEST(samples_csv, params.pep_dir, params.n_test)
    COMBINE_BFD_PROTEINS(BUILD_PEP_MANIFEST.out.manifest)
    MMSEQS_LINCLUST_BFD(COMBINE_BFD_PROTEINS.out.fasta.map { fasta -> tuple('bfd_proteins', fasta) })

    SPLIT_FASTA(
        MMSEQS_LINCLUST_BFD.out.representatives.map { prefix, fasta -> fasta },
        params.diamond_chunk_size,
    )
    chunks_ch = SPLIT_FASTA.out.chunks.flatten().map { f -> tuple(f.baseName, f) }

    DIAMOND_NOVELTY_SCREEN(chunks_ch)
    // .first() turns the single-emission dmnd channel into a value channel
    // so it broadcasts to every chunk rather than being consumed once and
    // starving DIAMOND_SEARCH_UNIPROT after the first chunk.
    DIAMOND_SEARCH_UNIPROT(chunks_ch, DIAMOND_MAKEDB.out.dmnd.first())

    BIN_NOVELTY_HITS(
        MMSEQS_LINCLUST_BFD.out.representatives.map { prefix, fasta -> fasta },
        DIAMOND_NOVELTY_SCREEN.out.map { chunk_id, tsv -> tsv }.collect(),
        DIAMOND_SEARCH_UNIPROT.out.map { chunk_id, tsv -> tsv }.collect(),
    )

    // Accumulation-curve analysis (design doc: BUILD_CLUSTER_INCIDENCE_MATRIX
    // + ACCUMULATION_CURVE). genome_classification_tsv is Stage 1's
    // (genome_classify.nf) output, a separate `nextflow run` invocation --
    // read here by path (params.genome_classification_tsv), the same
    // already-completed-output-by-path convention this file already uses for
    // params.uniprot_fungi_nr_fasta above, NOT a live inter-process channel
    // from a Stage-1 process (CLASSIFY_GENOMES's default, unnamed output
    // channel lives in genome_classify.nf's own workflow{}, a different
    // pipeline entirely).
    BUILD_CLUSTER_INCIDENCE_MATRIX(
        MMSEQS_LINCLUST_BFD.out.cluster_membership.map { prefix, tsv -> tsv },
        COMBINE_BFD_PROTEINS.out.provenance,
        BIN_NOVELTY_HITS.out.bins,
        samples_csv,
        Channel.fromPath(params.bfd_duckdb_busco_parquet, checkIfExists: true),
        Channel.fromPath(params.bfd_duckdb_asm_stats_parquet, checkIfExists: true),
        Channel.fromPath(params.genome_classification_tsv, checkIfExists: true),
    )

    ACCUMULATION_CURVE(
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.incidence_matrix,
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.external_status,
        BUILD_CLUSTER_INCIDENCE_MATRIX.out.genome_metadata,
        params.bfd_version_tag,
    )

    // Comparison-only branch (2026-09-08): dereplicate the same combined
    // BFD protein set a second way, then compare -- not fed into
    // SPLIT_FASTA/search/binning above.
    DIAMOND_CLUSTER(COMBINE_BFD_PROTEINS.out.fasta.map { fasta -> tuple('bfd_proteins', fasta) })
    EXTRACT_CLUSTER_REPS(DIAMOND_CLUSTER.out.clusters_and_fasta)

    COMPARE_CLUSTERINGS(
        MMSEQS_LINCLUST_BFD.out.cluster_membership.map { prefix, tsv -> tsv },
        EXTRACT_CLUSTER_REPS.out.cluster_membership.map { prefix, tsv -> tsv },
    )
}
