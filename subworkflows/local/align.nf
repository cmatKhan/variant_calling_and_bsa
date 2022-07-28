//
// Align reads to a reference genome
// note that this can be parameterized -- could put $param.aligner
// in the include ... from ... path below
//

include { BWAMEM2_ALIGNER                                 } from "${projectDir}/subworkflows/nf-core/align/bwamem2/main"
include { SAMTOOLS_SORT_INDEX_STATS                       } from "${projectDir}/subworkflows/nf-core/samtools/bam/sort_index_stats"
include { PICARD_ADDORREPLACEREADGROUPS as ADD_READ_GROUP } from "${projectDir}/modules/nf-core/modules/picard/addorreplacereadgroups/main"

workflow ALIGN {
    take:
    reads    // channel: [ val(meta), [ reads ] ]
    fasta   // channel: file(ref_genome)
    fasta_fai //channel: file(ref_genome) index
    intervals_bed_combined        // channel: [optional]  intervals_bed

    main:

    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()
    ch_bam      = Channel.empty()

    if(params.aligner == 'bwamem2') {
        BWAMEM2_ALIGNER (
            reads,
            fasta
        )
        ch_bam      = ch_bam.mix(BWAMEM2_ALIGNER.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_ALIGNER.out.versions)
    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }

    ADD_READ_GROUP(
        ch_bam
    )
    ch_versions = ch_versions.mix(ADD_READ_GROUP.out.versions)

    SAMTOOLS_SORT_INDEX_STATS(
        ADD_READ_GROUP.out.bam
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.mark_dups_report)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.stats)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.flagstat)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.idxstats)
    // Gather used softwares versions
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_STATS.out.versions)


    emit:
    bam_bai   = SAMTOOLS_SORT_INDEX_STATS.out.bam_index // channel: [ val(meta), [ bam ], [ bai ] ]
    report    = ch_reports                              // qc reports
    versions  = ch_versions                             // channel: [ versions.yml ]
}
