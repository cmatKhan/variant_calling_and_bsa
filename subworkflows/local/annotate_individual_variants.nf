//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF           } from "${projectDir}/modules/local/snpeff/main.nf"
include { TABIX_BGZIPTABIX } from "${projectDir}/modules/nf-core/modules/tabix/bgziptabix/main"

workflow ANNOTATE{
    take:
    vcf          // channel: [ val(meta), vcf ]

    main:
    ch_versions = Channel.empty()

    SNPEFF(vcf)
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = SNPEFF.out.report           //    path: *.csv
    versions = ch_versions                 //    path: versions.yml
}
