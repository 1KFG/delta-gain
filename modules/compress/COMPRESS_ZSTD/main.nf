// zstd-compress a large text output (FASTA/TSV) for publishing.
//
// Why a separate process, and why the plain file still exists in work/:
// mmseqs2, DIAMOND and the downstream python scripts in this pipeline read
// plain FASTA/TSV. None of the cached containers (python_3.12.sif, diamond,
// seqkit, mmseqs2) ships a zstd binary -- checked 2026-10-24 -- so producer
// processes keep emitting plain files for their consumers, and this process
// writes the compressed copy that is published to results/.
//
// Runs on the host (no container) with the shared-install zstd binary at
// params.zstd_bin. If that binary is missing, it falls back to pigz, then
// gzip. The output is integrity-tested before the task succeeds.

process COMPRESS_ZSTD {
    tag   { f.name }
    label 'compress_zstd'
    publishDir "${params.outdir}", mode: 'copy'

    input:
        path(f)

    output:
        path("${f.name}.{zst,gz}"), emit: compressed

    script:
    // zstd if params.zstd_bin is executable on the node; otherwise pigz;
    // otherwise gzip. The output suffix (.zst or .gz) shows which was used.
    """
    if [ -x "${params.zstd_bin}" ]; then
        ${params.zstd_bin} -T${task.cpus} -${params.zstd_level} -q -o ${f.name}.zst ${f}
        ${params.zstd_bin} -t -q ${f.name}.zst
    elif command -v pigz >/dev/null 2>&1; then
        pigz -p ${task.cpus} -6 -c ${f} > ${f.name}.gz
        pigz -t ${f.name}.gz
    else
        gzip -6 -c ${f} > ${f.name}.gz
        gzip -t ${f.name}.gz
    fi
    """

    stub:
    """
    touch ${f.name}.zst
    """
}
