
// include { MANTA_SOMATIC        } from "${projectDir}/modules/nf-core/modules/manta/somatic/main"
// include { CONTROLFREEC_SOMATIC          } from "${projectDir}/nf-core/variantcalling/controlfreec/somatic/main"
include { FREEBAYES as FREEBAYES_SINGLE } from "${projectDir}/modules/nf-core/modules/freebayes/main"
include { TIDDIT_SV                     } from "${projectDir}/modules/nf-core/modules/tiddit/sv/main.nf"

workflow CALL_INDIVIDUAL_VARIANTS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    fasta
    fasta_fai
    bwa_index // empty list passed for tiddit
    intervals                     // channel: [mandatory] intervals/target regions
    intervals_bed_gz_tbi          // channel: [mandatory] intervals/target regions index zipped and indexed
    intervals_bed_combined        // channel: [mandatory] intervals/target regions in one file unzipped

    main:

    ch_versions = Channel.empty()

    // MANTA_SOMATIC(
    //     bam_bai,
    //     fasta,
    //     fasta_fai
    //     )

    // manta_vcf                            = RUN_MANTA_SOMATIC.out.manta_vcf
    // manta_candidate_small_indels_vcf     = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf
    // manta_candidate_small_indels_vcf_tbi = RUN_MANTA_SOMATIC.out.manta_candidate_small_indels_vcf_tbi
    // ch_versions                          = ch_versions.mix(MANTA_SOMATIC.out.versions)

    // CONTROLFREEC_SOMATIC(
    //     controlfreec_input,
    //     fasta,
    //     fasta_fai,
    //     dbsnp,
    //     dbsnp_tbi,
    //     chr_files,
    //     mappability,
    //     intervals_bed_combined
    // )

    // ch_versions = ch_versions.mix(CONTROLFREEC_SOMATIC.out.versions)

    freebayes_input = bam_bai
            .map{ meta, bam, bai, intervals ->
                [meta, bam, bai, [], [], intervals]
            }

    FREEBAYES_SINGLE(
        freebayes_input,
        fasta,
        fasta_fai
    )
    freebayes_vcf = FREEBAYES_SINGLE.out.freebayes_vcf
    ch_versions   = ch_versions.mix(FREEBAYES_SINGLE.out.versions)

    TIDDIT_SV(
        bam_bai,
        fasta,
        bwa_index
    )

    tiddit_ploidy = TIDDIT_SV.out.ploidy

    tiddit_vcf = TIDDIT_SV.out.tiddit_vcf
    ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)


    emit:
    freebayes_vcf           // channel: [val(meta), path(vcf) ]
    tiddit_vcf              // channel: [val(meta), path(vcf) ]
    tiddit_ploidy = TIDDIT_SV.out.ploidy
    versions  = ch_versions // channel: [ versions.yml ]
}
