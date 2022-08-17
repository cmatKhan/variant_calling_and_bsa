
include { TIDDIT_SV } from "${projectDir}/modules/nf-core/modules/tiddit/sv/main.nf"
include { VCFTOOLS        } from "${projectDir}/modules/nf-core/modules/vcftools/main"

workflow CALL_STRUCTURAL_VARIANTS {
    take:
    bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    dusted_bed
    fasta
    bwa_index

    main:

    ch_versions = Channel.empty()


    TIDDIT_SV(
        bam_bai,
        fasta,
        bwa_index
    )
    ch_versions = ch_versions.mix(TIDDIT_SV.out.versions)

    VCFTOOLS(
        TIDDIT_SV.out.vcf,
        dusted_bed,
        []
    )
    ch_versions   = ch_versions.mix(VCFTOOLS.out.versions)


    emit:
    vcf       = VCFTOOLS.out.vcf       // channel: [val(meta), path(vcf) ]
    ploidy    = TIDDIT_SV.out.ploidy
    versions  = ch_versions       // channel: [ versions.yml ]
}
