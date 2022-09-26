
// include { MANTA_SOMATIC        } from "${projectDir}/modules/nf-core/modules/manta/somatic/main"
// include { CONTROLFREEC_SOMATIC          } from "${projectDir}/nf-core/variantcalling/controlfreec/somatic/main"
include { FREEBAYES as FREEBAYES_INDIVIDUAL } from "${projectDir}/modules/nf-core/modules/freebayes/main"
include { CNVPYTOR_IMPORTREADDEPTH          } from "${projectDir}/modules/nf-core/modules/cnvpytor/importreaddepth/main"
include { CNVPYTOR_HISTOGRAM                } from "${projectDir}/modules/nf-core/modules/cnvpytor/histogram/main"
include { CNVPYTOR_PARTITION                } from "${projectDir}/modules/nf-core/modules/cnvpytor/partition/main"
include { CNVPYTOR_CALLCNVS                 } from "${projectDir}/modules/nf-core/modules/cnvpytor/callcnvs/main"
include { CNVPYTOR_VIEW                     } from "${projectDir}/modules/nf-core/modules/cnvpytor/view/main"
include { VCFTOOLS                          } from "${projectDir}/modules/nf-core/modules/vcftools/main"

workflow CALL_INDIVIDUAL_VARIANTS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    dusted_bed
    fasta
    fasta_fai  // channel: [val(meta), path(fai)]
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped
    cnvpytor_genome_conf
    cnvpytor_genome_gc_ch
    cnv_histogram_bin_size
    cnv_output_format

    main:

    ch_versions = Channel.empty()
    freebayes_vcf = Channel.empty()

    fasta_fai = fasta_fai.map(it -> it[1])

    bam_bai
        .map{ meta, bam, bai ->
            [meta, bam, bai, [], []]}
        .combine(intervals_bed_combined)
        .set{ freebayes_input }

    FREEBAYES_INDIVIDUAL(
        freebayes_input,
        fasta,
        fasta_fai,
        [],
        [],
        []
    )
    ch_versions   = ch_versions.mix(FREEBAYES_INDIVIDUAL.out.versions)

    // create channel to collect vcf files from freebayes and cnvnator
    FREEBAYES_INDIVIDUAL.out.vcf.set{ ch_vcf }

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
        cnv_output_format.map(it -> it[0]),
        cnvpytor_genome_conf,
        cnvpytor_genome_gc_ch
    )
    ch_versions = ch_versions.mix(CNVPYTOR_VIEW.out.versions)

    ch_vcf = ch_vcf.mix(CNVPYTOR_VIEW.out.vcf)

    VCFTOOLS(
        ch_vcf,
        dusted_bed,
        []
    )
    ch_versions   = ch_versions.mix(VCFTOOLS.out.versions)



    emit:
    freebayes_vcf = VCFTOOLS.out.vcf // channel: [val(meta), path(vcf) ]
    versions  = ch_versions          // channel: [ versions.yml ]
}
