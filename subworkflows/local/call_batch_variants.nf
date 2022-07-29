
include { FREEBAYES_BATCH } from "${projectDir}/modules/local/freebayes_batch/main"

workflow CALL_BATCH_VARIANTS {
    take:
    freebayes_input  // channel: [ val(meta), path(bam), path(bai), path(target_bed) ]
    fasta
    fasta_fai  // channel: [val(meta), path(fai)]

    main:

    ch_versions = Channel.empty()
    freebayes_vcf = Channel.empty()
    tiddit_ploidy = Channel.empty()
    tiddit_vcf = Channel.empty()
    tiddit_ploidy = Channel.empty()

    fasta_fai = fasta_fai.map(it -> it[1])

    FREEBAYES_BATCH(
        freebayes_input,
        fasta,
        fasta_fai
    )
    freebayes_vcf = FREEBAYES_BATCH.out.vcf
    ch_versions   = ch_versions.mix(FREEBAYES_BATCH.out.versions)

    emit:
    freebayes_vcf = FREEBAYES_BATCH.out.vcf           // channel: [val(meta), path(vcf) ]
    tiddit_vcf              // channel: [val(meta), path(vcf) ]
    tiddit_ploidy
    versions  = ch_versions // channel: [ versions.yml ]
}
