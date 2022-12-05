
include { CNVPYTOR_IMPORTREADDEPTH          } from "${projectDir}/modules/nf-core/modules/cnvpytor/importreaddepth/main"
include { CNVPYTOR_HISTOGRAM                } from "${projectDir}/modules/nf-core/modules/cnvpytor/histogram/main"
include { CNVPYTOR_PARTITION                } from "${projectDir}/modules/nf-core/modules/cnvpytor/partition/main"
include { CNVPYTOR_CALLCNVS                 } from "${projectDir}/modules/nf-core/modules/cnvpytor/callcnvs/main"
include { CNVPYTOR_VIEW                     } from "${projectDir}/modules/nf-core/modules/cnvpytor/view/main"

workflow CNVPYTOR_COMPLETE {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    cnvpytor_genome_conf
    cnvpytor_genome_gc_ch
    cnv_histogram_bin_size
    cnv_output_format

    CNVPYTOR_IMPORTREADDEPTH(
        bam_bai,
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_IMPORTREADDEPTH.out.versions)

    cnv_histogram_bin_size.map{it -> it[0]}.set{ bin_size }

    CNVPYTOR_HISTOGRAM(
        CNVPYTOR_IMPORTREADDEPTH.out.pytor,
        bin_size,
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_HISTOGRAM.out.versions)

    CNVPYTOR_PARTITION(
        CNVPYTOR_HISTOGRAM.out.pytor,
        bin_size,
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_PARTITION.out.versions)

    CNVPYTOR_CALLCNVS(
        CNVPYTOR_PARTITION.out.pytor,
        bin_size,
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_CALLCNVS.out.versions)

    CNVPYTOR_VIEW(
        CNVPYTOR_PARTITION.out.pytor,
        bin_size,
        cnv_output_format,
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_VIEW.out.versions)
}
