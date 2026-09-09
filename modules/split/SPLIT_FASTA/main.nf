// Chunk the dereplicated BFD-proteins representative FASTA into fixed-size
// pieces so DIAMOND_NOVELTY_SCREEN can fan out across SLURM tasks instead of
// running as one giant serial search -- same reasoning as genome_classify's
// per-component SKANI_TRIANGLE fan-out, just chunked by sequence count here
// rather than by mash component.

process SPLIT_FASTA {
    label 'split_fasta'

    input:
        path(fasta)
        val(chunk_size)

    output:
        path("chunks/*.fasta"), emit: chunks

    script:
    """
    mkdir -p chunks
    seqkit split2 -s ${chunk_size} -O chunks -o chunk ${fasta}
    """

    stub:
    """
    mkdir -p chunks
    printf ">STUB_Q1 stub\\nMSTUBSEQ\\n" > chunks/chunk.part_001.fasta
    """
}
