
// include { MANTA_SOMATIC        } from "${projectDir}/modules/nf-core/modules/manta/somatic/main"
// include { CONTROLFREEC_SOMATIC          } from "${projectDir}/nf-core/variantcalling/controlfreec/somatic/main"
include { FREEBAYES as FREEBAYES_INDIVIDUAL } from "${projectDir}/modules/nf-core/modules/freebayes/main"
include { VCFTOOLS                          } from "${projectDir}/modules/nf-core/modules/vcftools/main"

workflow CALL_INDIVIDUAL_VARIANTS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    dusted_bed
    fasta
    fasta_fai  // channel: [val(meta), path(fai)]
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped

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

    VCFTOOLS(
        FREEBAYES_INDIVIDUAL.out.vcf,
        dusted_bed,
        []
    )
    ch_versions   = ch_versions.mix(VCFTOOLS.out.versions)



    emit:
    freebayes_vcf = VCFTOOLS.out.vcf // channel: [val(meta), path(vcf) ]
    versions  = ch_versions          // channel: [ versions.yml ]
}
