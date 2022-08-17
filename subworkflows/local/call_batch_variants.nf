
include { FREEBAYES_BATCH } from "${projectDir}/modules/local/freebayes_batch/main"
include { GATK4_MERGEVCFS } from "${projectDir}/modules/nf-core/modules/gatk4/mergevcfs/main"
include { VCFTOOLS        } from "${projectDir}/modules/nf-core/modules/vcftools/main"

workflow CALL_BATCH_VARIANTS {
    take:
    freebayes_input  // channel: [ val(meta), path(bam), path(bai), path(target_bed) ]
    dusted_bed
    fasta
    fasta_fai  // channel: [val(meta), path(fai)]
    sequence_dict

    main:

    ch_versions = Channel.empty()

    fasta_fai = fasta_fai.map(it -> it[1])

    FREEBAYES_BATCH(
        freebayes_input,
        fasta,
        fasta_fai
    )
    ch_versions   = ch_versions.mix(FREEBAYES_BATCH.out.versions)

    freebayes_combined =
        FREEBAYES_BATCH.out.vcf
        .map{ meta, vcf ->
            [meta.id, meta.aligner, vcf]}
        .groupTuple(by: [0,1])
        .map{ group, aligner, vcf_list ->
            [[id:group, aligner: aligner], vcf_list]}
        .set{ ch_collected_vcf }

    GATK4_MERGEVCFS( ch_collected_vcf, sequence_dict )
    ch_versions   = ch_versions.mix(GATK4_MERGEVCFS.out.versions)

    VCFTOOLS(
        GATK4_MERGEVCFS.out.vcf,
        dusted_bed,
        []
    )
    ch_versions   = ch_versions.mix(VCFTOOLS.out.versions)


    emit:
    freebayes_vcf = VCFTOOLS.out.vcf // channel: [val(meta), path(vcf) ]
    versions  = ch_versions          // channel: [ versions.yml ]
}
