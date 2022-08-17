//
// Align reads to a reference genome
// note that this can be parameterized -- could put $param.aligner
// in the include ... from ... path below
//

include { BWAMEM2_ALIGNER               } from "${projectDir}/subworkflows/nf-core/align/bwamem2/main"
include { NGM_YAHA                      } from "${projectDir}/subworkflows/nf-core/align/nextgenmap_yaha/main"
include { SAMTOOLS_SORT_INDEX_STATS     } from "${projectDir}/subworkflows/nf-core/samtools/bam/sort_index_stats"
include { PICARD_ADDORREPLACEREADGROUPS } from "${projectDir}/modules/nf-core/modules/picard/addorreplacereadgroups/main"

workflow ALIGN {
    take:
    aligners
    reads    // channel: [ val(meta), [ reads ] ]
    fasta    // channel: path(fasta)
    fai
    bwamem2_index
    yaha_index
    yaha_nib2

    main:

    ch_versions    = Channel.empty()
    ch_reports     = Channel.empty()
    ch_bam         = Channel.empty()

    def aligners_unrecognized = true

    if(aligners && aligners.split(',').contains('bwamem2')) {
        aligners_unrecognized = false

        reads.map{meta, reads ->
            def meta_tmp = augment_meta(meta, "bwamem2")
                [meta_tmp,reads]}
        .set{bwamem2_input}

        BWAMEM2_ALIGNER (
            bwamem2_input,
            bwamem2_index
        )
        ch_bam      = ch_bam.mix(BWAMEM2_ALIGNER.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_ALIGNER.out.versions)
    }

    if(aligners && aligners.split(',').contains('ngm_yaha')){
        aligners_unrecognized = false

        reads.map{meta, reads ->
            def meta_tmp = augment_meta(meta, "ngm_yaha")
                [meta_tmp,reads]}
        .set{ngm_yaha_input}

        NGM_YAHA (
            ngm_yaha_input,
            fasta,
            fai,
            yaha_index,
            yaha_nib2
        )
        ch_bam      = ch_bam.mix(NGM_YAHA.out.bam)
        ch_versions = ch_versions.mix(NGM_YAHA.out.versions)
    }

    // if no aligner block is run, then exit with error
    if(aligners_unrecognized){
        exit 1, "Either no aligners were specified, or no specified aligner " +
        "is recognized. Aligners in list: ${aligners.split(',')}."
    }

    PICARD_ADDORREPLACEREADGROUPS(
        ch_bam
    )
    ch_versions = ch_versions.mix(PICARD_ADDORREPLACEREADGROUPS.out.versions)

    SAMTOOLS_SORT_INDEX_STATS(
        PICARD_ADDORREPLACEREADGROUPS.out.bam
    )

    // Gather QC reports
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.mark_dups_report)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.stats)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.flagstat)
    ch_reports  = ch_reports.mix(SAMTOOLS_SORT_INDEX_STATS.out.idxstats)
    // Gather used softwares versions
    ch_versions = ch_versions.mix(SAMTOOLS_SORT_INDEX_STATS.out.versions)

    emit:
    // channel: [ val(meta), [ bam ], [ bai ] ]
    // meta includes aligner: <aligner>, eg aligner: bwamem2
    bam_bai   = SAMTOOLS_SORT_INDEX_STATS.out.bam_index
    report    = ch_reports                              // qc reports
    versions  = ch_versions                             // channel: [ versions.yml ]
}

def augment_meta(Map meta, aligner) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta["aligner"] = aligner

    return new_meta
}
