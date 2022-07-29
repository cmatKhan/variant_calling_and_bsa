//
// Index the reference genome and align reads
//

include { NEXTGENMAP           } from "${projectDir}/modules/local/nextgenmap/main"
include { NEXTGENMAP_PARTITION } from "${projectDir}/subworkflows/nf-core/align/nextgenmap_yaha/nextgenmap_partition.nf"
include { YAHA_ALIGN           } from "${projectDir}/modules/local/yaha/align/main"

workflow NGM_YAHA {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    fasta // file(fasta)
    fai
    yaha_index
    yaha_nib2

    main:

    ch_versions = Channel.empty()

    NEXTGENMAP( reads, fasta )
    ch_versions = ch_versions.mix(NEXTGENMAP.out.versions.first())

    NEXTGENMAP_PARTITION( NEXTGENMAP.out.bam, fasta )
    ch_versions = ch_versions.mix(NEXTGENMAP_PARTITION.out.versions.first())

    fai.map{meta,fai -> return fai}
        .set{ samtools_merge_fai }

    YAHA_ALIGN(
        NEXTGENMAP_PARTITION.out.unaligned_reads,
        NEXTGENMAP_PARTITION.out.aligned,
        yaha_nib2,
        yaha_index,
        fasta,
        samtools_merge_fai
    )
    ch_versions = ch_versions.mix(YAHA_ALIGN.out.versions.first())

    emit:
    bam            = YAHA_ALIGN.out.ngm_yaha_bam  // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions // channel: [ versions.yml ]
}
