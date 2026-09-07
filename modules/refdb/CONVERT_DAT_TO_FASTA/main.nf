// Convert a UniProt .dat.gz flatfile to FASTA + an InterPro/Pfam/GO
// cross-reference table (bin/dat_to_fasta_and_annotations.py).
//
// Requires biopython/1.87, NOT the default python3's biopython (1.79):
// live-tested 2026-09-06 -- Bio.SwissProt.parse under 1.79 throws
// `IndexError: list index out of range` in its FT (feature table) line
// parser on current (2026-release) UniProt flatfile syntax; 1.87 parses the
// identical real record correctly. Pin the module explicitly, inline in the
// script (per this cluster's beforeScript-doesn't-reach-script: gotcha) --
// don't rely on whatever biopython the login shell's default python3 has.
//
// module load (not apptainer) deliberately here: this runs twice, serially,
// no concurrency contention -- unlike skani/mash's squashfuse-under-load
// failure mode (see profile_ANI.config), there's no concurrency problem here
// to solve by containerizing.

process CONVERT_DAT_TO_FASTA {
    tag   { division }
    label 'convert_dat'

    input:
        tuple val(division), path(dat_gz)

    output:
        tuple val(division), path("${division}_fungi.fasta"), path("${division}_fungi.xrefs.tsv")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load biopython/1.87
    ${projectDir}/bin/dat_to_fasta_and_annotations.py \\
        --dat ${dat_gz} \\
        --prefix ${division}_fungi
    """

    stub:
    """
    printf ">STUB1 stub\\nMSTUBSEQ\\n" > ${division}_fungi.fasta
    printf "accession\\tdb\\tid\\textra1\\textra2\\n" > ${division}_fungi.xrefs.tsv
    """
}
