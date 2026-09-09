// Dereplicate a protein set via MMseqs2 linclust -- the linear-time
// clustering algorithm, the only one of MMseqs2's clustering modes that
// scales past low tens of millions of sequences (the design doc's Step 2
// already flagged cascaded `mmseqs cluster` as not scaling to
// hundreds-of-millions-of-sequence inputs; the same reasoning applies here).
//
// Generalized 2026-09-08 to take an explicit `prefix` (originally hardcoded
// to "uniprot_fungi_*"): reused as-is for BFD's own ~176M-protein set
// (Stage 3 BFD-side dereplication), not just the UniProt-fungi reference
// set this was first written for -- same tool, same linclust settings
// already proven at UniProt-fungi scale (see HANDOFF.md item 3).
//
// Output: a representative-sequence FASTA (one sequence per cluster) plus a
// cluster membership TSV (representative accession <TAB> member accession),
// so cluster-to-original-accession provenance is preserved for joining back
// to per-source annotation tables downstream.

process MMSEQS_LINCLUST {
    tag   { prefix }
    label 'mmseqs_linclust'
    // Missing until 2026-09-07: the real run completed successfully but its
    // output only existed in work/ (never copied out) until this was added
    // -- caught after the fact, output recovered manually from the work dir
    // for that run. Don't repeat for future runs.
    publishDir "${params.outdir}", mode: 'copy'

    input:
        tuple val(prefix), path(combined_fasta)

    output:
        tuple val(prefix), path("${prefix}_nr.fasta"),      emit: representatives
        tuple val(prefix), path("${prefix}_clusters.tsv"),  emit: cluster_membership

    script:
    """
    mkdir -p mmseqs_tmp
    mmseqs createdb ${combined_fasta} DB

    mmseqs linclust DB DB_clu mmseqs_tmp \\
        --min-seq-id ${params.mmseqs_min_seq_id} \\
        -c ${params.mmseqs_coverage} \\
        --cov-mode ${params.mmseqs_cov_mode} \\
        --threads ${task.cpus}

    mmseqs createsubdb DB_clu DB DB_rep
    mmseqs convert2fasta DB_rep ${prefix}_nr.fasta

    mmseqs createtsv DB DB DB_clu ${prefix}_clusters.tsv
    """

    stub:
    """
    printf ">STUB_REP1 stub\\nMSTUBSEQ\\n" > ${prefix}_nr.fasta
    printf "STUB_REP1\\tSTUB_REP1\\n" > ${prefix}_clusters.tsv
    """
}
