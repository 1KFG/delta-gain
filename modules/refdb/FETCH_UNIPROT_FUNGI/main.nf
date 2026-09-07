// Stage 3 reference-database staging: pull UniProtKB fungi (taxid 4751) as
// FASTA via the REST stream endpoint. This is the "quick and dirty" v0
// baseline the design doc converged on (see DeltaGain/README.md Stage 3 --
// nr_cluster_seq.dmnd covers the broad novelty screen; this pulls the
// specifically-fungal, UniProt-curated baseline for the taxonomic-breadth
// axis and as a smaller/cleaner complement to the nr_cluster screen).
//
// Live-tested 2026-09-06 against the real API (small record counts only --
// see README.md for what was and wasn't pulled). Two real findings from that
// test, both load-bearing for how this process is written:
//
//   1. The /stream endpoint has a hard SERVER-SIDE result-count cap. A
//      `reviewed:true` (Swiss-Prot-only) fungi query streams fine. The
//      unrestricted fungi query (Swiss-Prot + TrEMBL, taxid 4751 with no
//      reviewed filter) is rejected outright with HTTP 403 "Too many results
//      to retrieve. Please refine your query or consider fetching batch by
//      batch" -- BEFORE any data streams. This is not something client-side
//      truncation (piping through awk/head) can work around: the server
//      refuses based on the query's total result count, not on how much the
//      client asks to read.
//   2. There is no server-side small-N knob either way -- passing `size=` to
//      /stream returns the same 403, even for size=1.
//
// Consequence: 'reviewed' mode (this file's only currently-implemented mode)
// is real and safe to run at any record count. A full Swiss-Prot+TrEMBL
// fungi pull needs PAGINATION (UniProt's /uniprotkb/search endpoint, cursor-
// based via the response `Link` header, looping in pages of ~500) -- a
// materially different implementation, not built here. Deliberately NOT
// built in this session (only a small connectivity/content test was wanted,
// not the real full pull) -- see mode: 'all_unpaginated' below, which exists
// only to fail fast with this explanation rather than silently attempt (and
// waste retries on) a call that is guaranteed to 403.

process FETCH_UNIPROT_FUNGI {
    tag   { mode }
    label 'fetch_uniprot'
    publishDir "${params.refdb_outdir}", mode: 'copy'

    input:
        val(mode)   // 'reviewed' (Swiss-Prot only, works) | 'all_unpaginated' (fails fast, see header comment)

    output:
        path("uniprot_fungi_${mode}.fasta"), emit: fasta

    script:
    if (mode == 'all_unpaginated') {
        """
        echo "FETCH_UNIPROT_FUNGI mode=all_unpaginated: not implemented -- the" >&2
        echo "unrestricted Swiss-Prot+TrEMBL fungi query exceeds UniProt's" >&2
        echo "/stream endpoint result-count cap (confirmed live 2026-09-06:" >&2
        echo "HTTP 403 'Too many results to retrieve... fetch batch by batch')." >&2
        echo "Needs cursor-based pagination via /uniprotkb/search instead --" >&2
        echo "see this module's header comment. Not built yet." >&2
        exit 1
        """
    } else {
        def query = "taxonomy_id:${params.fungi_taxid}+AND+reviewed:true"
        def truncate = (params.uniprot_test_records as int) > 0
            ? "awk '/^>/{c++} c>${params.uniprot_test_records as int}{exit} {print}'"
            : "cat"
        """
        curl -fsS --retry 3 --retry-delay 5 \\
            "https://rest.uniprot.org/uniprotkb/stream?query=${query}&format=fasta" \\
            | ${truncate} > uniprot_fungi_${mode}.fasta

        n=\$(grep -c '^>' uniprot_fungi_${mode}.fasta || true)
        echo "FETCH_UNIPROT_FUNGI: wrote \${n} records to uniprot_fungi_${mode}.fasta" >&2
        if [ "\${n:-0}" -eq 0 ]; then
            echo "FETCH_UNIPROT_FUNGI: zero records written -- treating as a failure" >&2
            exit 1
        fi
        """
    }

    stub:
    """
    printf ">sp|STUB1|STUB_FUNGI stub record\\nMSTUBSEQ\\n" > uniprot_fungi_${mode}.fasta
    """
}
