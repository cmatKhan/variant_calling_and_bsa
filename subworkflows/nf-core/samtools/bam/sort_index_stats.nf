//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
// COPIED FROM nf-co/rnaseq
//

include { SAMTOOLS_SORT         } from "${projectDir}/modules/nf-core/modules/samtools/sort/main"
include { SAMTOOLS_INDEX        } from "${projectDir}/modules/nf-core/modules/samtools/index/main"
include { SAMTOOLS_BAM_STATS    } from "${projectDir}/subworkflows/nf-core/samtools/bam/bam_stats"
include { PICARD_MARKDUPLICATES } from "${projectDir}/modules/nf-core/modules/picard/markduplicates/main"

workflow SAMTOOLS_SORT_INDEX_STATS {
    take:
    bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_SORT ( bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    PICARD_MARKDUPLICATES(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    PICARD_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_sorted_bam_index }

    SAMTOOLS_BAM_STATS ( ch_sorted_bam_index )
    ch_versions = ch_versions.mix(SAMTOOLS_BAM_STATS.out.versions)

    emit:
    bam_index  = ch_sorted_bam_index             // channel: [ val(meta), path(bam), path(bai) ]
    mark_dups_report = PICARD_MARKDUPLICATES.out.metrics // channel: mark dups report
    stats      = SAMTOOLS_BAM_STATS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat   = SAMTOOLS_BAM_STATS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats   = SAMTOOLS_BAM_STATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    versions   = ch_versions                     // channel: [ versions.yml ]
}
