// Bulk fetch of UniProt's pre-split taxonomic-division fungi files (FTP),
// NOT the REST /stream endpoint -- see FETCH_UNIPROT_FUNGI's header comment
// for why /stream can't serve the full (comprehensive) fungi set: it 403s
// outright on the unrestricted Swiss-Prot+TrEMBL query, before any data
// streams, regardless of client-side truncation.
//
// Live-verified 2026-09-06 (HEAD request, see conversation): both files exist
// at this path for the current release.
//   uniprot_sprot_fungi.dat.gz    ~60 MB   (Swiss-Prot, reviewed, curated)
//   uniprot_trembl_fungi.dat.gz   ~13 GB   (TrEMBL, unreviewed, automatic)
// Together these are the comprehensive UniProtKB fungi set (confirmed via
// REST total-count headers: 38,451 reviewed + 17,534,069 unreviewed =
// 17,572,520, matching the unrestricted taxonomy_id:4751 total exactly).
//
// -C - (curl resume) matters here: the trembl file is large enough that a
// dropped connection losing all progress would be a real cost, not a
// nuisance.

process FETCH_UNIPROT_FUNGI_FTP {
    tag   { division }
    label 'fetch_uniprot_ftp'
    publishDir "${params.uniprot_staging_dir}", mode: 'copy'

    input:
        tuple val(division), val(remote_name)   // division: 'sprot' | 'trembl'

    output:
        tuple val(division), path(remote_name), emit: dat_gz

    script:
    """
    curl -fsSL --retry 5 --retry-delay 30 -C - \\
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/${remote_name}" \\
        -o ${remote_name}

    sz=\$(stat -c%s ${remote_name} 2>/dev/null || stat -f%z ${remote_name})
    echo "FETCH_UNIPROT_FUNGI_FTP: ${remote_name} staged, \${sz} bytes" >&2
    if [ "\${sz:-0}" -lt 1000000 ]; then
        echo "FETCH_UNIPROT_FUNGI_FTP: suspiciously small (<1MB) -- treating as a failure" >&2
        exit 1
    fi
    """

    stub:
    """
    printf "stub placeholder, not a real gzip" > ${remote_name}
    """
}
