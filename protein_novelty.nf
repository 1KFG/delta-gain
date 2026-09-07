#!/usr/bin/env nextflow
/*
 * protein_novelty.nf — DeltaGain Stage 3 (in progress).
 *
 * Two independent pieces wired so far:
 *   1. FETCH_UNIPROT_FUNGI (REST /stream, reviewed-only) — small connectivity/
 *      content smoke test. Live-tested 2026-09-06.
 *   2. FETCH_UNIPROT_FUNGI_FTP -> CONVERT_DAT_TO_FASTA -> COMBINE_UNIPROT_FASTA
 *      -> MMSEQS_LINCLUST — the real, comprehensive UniProtKB-fungi
 *      (Swiss-Prot + TrEMBL, ~17.5M sequences) reference set, bulk-downloaded
 *      from FTP (not REST -- see FETCH_UNIPROT_FUNGI_FTP's header comment for
 *      why /stream can't serve this), converted to FASTA + InterPro/Pfam/GO
 *      cross-refs, combined, and dereplicated via MMseqs2 linclust.
 *
 * NOT YET RUN FOR REAL past small tests: the FTP path pulls a 13 GB file
 * (uniprot_trembl_fungi.dat.gz) and linclust-clusters ~17.5M sequences -- a
 * genuinely long-running, resource-heavy step. Stub-validated only until
 * explicitly run for real.
 *
 * Still not wired in: DIAMOND_NOVELTY_SCREEN (BFD proteins vs. nr_cluster_seq.dmnd,
 * module already written), and BFD-proteins-vs-this-UniProt-set search --
 * this workflow currently only builds the UniProt-fungi reference set itself.
 */

nextflow.enable.dsl = 2

include { FETCH_UNIPROT_FUNGI }     from './modules/refdb/FETCH_UNIPROT_FUNGI/main.nf'
include { FETCH_UNIPROT_FUNGI_FTP } from './modules/refdb/FETCH_UNIPROT_FUNGI_FTP/main.nf'
include { CONVERT_DAT_TO_FASTA }    from './modules/refdb/CONVERT_DAT_TO_FASTA/main.nf'
include { COMBINE_UNIPROT_FASTA }   from './modules/refdb/COMBINE_UNIPROT_FASTA/main.nf'
include { MMSEQS_LINCLUST }         from './modules/dedup/MMSEQS_LINCLUST/main.nf'

workflow {
    // Small connectivity/content smoke test (REST) -- cheap, safe at any scale.
    FETCH_UNIPROT_FUNGI(Channel.of('reviewed'))

    // The real comprehensive reference set (FTP bulk download).
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

    MMSEQS_LINCLUST(COMBINE_UNIPROT_FASTA.out.fasta)
}
