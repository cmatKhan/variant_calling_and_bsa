//
// Run SNPEFF to annotate VCF files
//

include { SNPEFF           } from "${projectDir}/modules/local/snpeff/main.nf"
include { TABIX_BGZIPTABIX } from "${projectDir}/modules/nf-core/modules/tabix/bgziptabix/main"
include { PICARD_SORTVCF   } from "${projectDir}/modules/nf-core/modules/picard/sortvcf/main"

workflow ANNOTATE{
    take:
    vcf          // channel: [ val(meta), path(vcf) ]
    fasta       // channel: path(fasta file)

    main:
    ch_versions = Channel.empty()
    ch_report = Channel.empty()
    snpeff_config = Channel.fromPath(params.snpeff_db_config).collect()
    snpeff_db = Channel.fromPath("${projectDir}/assets/snpeff_db").collect()

    PICARD_SORTVCF(vcf,fasta,[])
    SNPEFF(PICARD_SORTVCF.out.vcf, snpeff_config, snpeff_db)
    TABIX_BGZIPTABIX(SNPEFF.out.vcf)

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(PICARD_SORTVCF.out.versions.first())
    ch_versions = ch_versions.mix(SNPEFF.out.versions.first())
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    ch_report = ch_report.mix(SNPEFF.out.report)

    emit:
    vcf_tbi  = TABIX_BGZIPTABIX.out.gz_tbi // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    reports  = ch_report                   //    path: *.csv
    versions = ch_versions                 //    path: versions.yml
}
