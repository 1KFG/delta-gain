// ONE global all-vs-all mash comparison across every BFD genome in the run,
// with no prior taxonomic bucketing. This is the label-free step that lets a
// genuinely novel/mislabeled species be discovered by distance rather than by
// whatever species name it happened to be deposited under.
//
// BFD-internal only, deliberately: this stage answers "which BFD genomes are
// redundant strains of each other," which needs no external reference genome.
// (An earlier draft of this pipeline mixed UniProt/RefSeq source genomes into
// this same comparison — unnecessary: whether a species already has any
// annotation elsewhere is a taxonomy lookup, not a genome-ANI question, and
// the protein-level novelty screen in Stage 3 needs a protein reference
// database, not genome assemblies, and doesn't care which assembly a known
// protein came from. See DeltaGain/README.md.)
//
// mash dist is an approximate (MinHash) distance, cheap enough to run at full
// 22k-genome scale in one shot; skani triangle refinement only happens
// afterwards, per connected component (see FIND_COMPONENTS + SKANI_TRIANGLE),
// so the expensive exact-ANI step never has to be all-vs-all.

include { mashCpusFor; mashMemoryFor; mashTimeFor } from '../../../lib/utils.nf'

process MASH_PREFILTER_GLOBAL {
    tag    "all n=${sketches.size()}"
    label  'mash_prefilter'

    cpus   { mashCpusFor(sketches.size()) }
    memory { mashMemoryFor(sketches.size(), task.attempt) }
    time   { mashTimeFor(sketches.size(), task.attempt) }

    input:
        path(sketches)

    output:
        path("all_genomes.mash_dist.tsv")

    script:
    """
    ls *.msh > sketch_list.txt
    mash paste combined -l sketch_list.txt
    mash dist -p ${task.cpus} combined.msh combined.msh > all_genomes.mash_dist.tsv
    """

    stub:
    """
    touch all_genomes.mash_dist.tsv
    """
}
