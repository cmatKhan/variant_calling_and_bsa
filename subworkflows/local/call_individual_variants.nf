
// include { MANTA_SOMATIC        } from "${projectDir}/modules/nf-core/modules/manta/somatic/main"
// include { CONTROLFREEC_SOMATIC          } from "${projectDir}/nf-core/variantcalling/controlfreec/somatic/main"
include { FREEBAYES as FREEBAYES_INDIVIDUAL } from "${projectDir}/modules/nf-core/modules/freebayes/main"
include { TIDDIT_SV                         } from "${projectDir}/modules/nf-core/modules/tiddit/sv/main.nf"

workflow CALL_INDIVIDUAL_VARIANTS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    fasta
    fasta_fai  // channel: [val(meta), path(fai)]
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped

    main:

    ch_versions = Channel.empty()
    freebayes_vcf = Channel.empty()
    tiddit_ploidy = Channel.empty()
    tiddit_vcf = Channel.empty()
    tiddit_ploidy = Channel.empty()

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
    freebayes_vcf = FREEBAYES_INDIVIDUAL.out.vcf
    ch_versions   = ch_versions.mix(FREEBAYES_INDIVIDUAL.out.versions)

    // TIDDIT_SV(
    //     bam_bai,
    //     fasta,
    //     bwa_index
    // )

    // tiddit_ploidy = TIDDIT_SV.out.ploidy

    // tiddit_vcf = TIDDIT_SV.out.vcf
    // ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)


    emit:
    freebayes_vcf = FREEBAYES_INDIVIDUAL.out.vcf           // channel: [val(meta), path(vcf) ]
    tiddit_vcf              // channel: [val(meta), path(vcf) ]
    tiddit_ploidy
    versions  = ch_versions // channel: [ versions.yml ]
}
