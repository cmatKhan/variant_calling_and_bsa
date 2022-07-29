//
// parition ngm alignments to primary alignments, and unmapped reads.
// filter out unmapped reads which do not fail vendor QC (0x512) from NGM
// alignments. Merge unmapped R1 and R2 with Pear, and cat to create a 'pseudo'
// single end fastq file for YAHA
//

include { SAMTOOLS_SORT                   } from "${projectDir}/modules/nf-core/modules/samtools/sort/main"
include { PICARD_MARKDUPLICATES           } from "${projectDir}/modules/nf-core/modules/picard/markduplicates/main.nf"
include { SAMTOOLS_INDEX                  } from "${projectDir}/modules/nf-core/modules/samtools/index/main"
include { SAMTOOLS_VIEW as EXTRACT_MAPPED } from "${projectDir}/modules/nf-core/modules/samtools/view/main"
include { SAMTOOLS_FASTQ                  } from "${projectDir}/modules/nf-core/modules/samtools/fastq/main"
include { EXTRACT_UMAP                    } from "${projectDir}/modules/local/extract_umap.nf"
include { PEAR                            } from "${projectDir}/modules/nf-core/modules/pear/main"
include { CAT_FASTQ                       } from "${projectDir}/modules/nf-core/modules/cat/fastq/main"


workflow NEXTGENMAP_PARTITION {
    take:
    bam // channel: [ val(meta), [ bam ] ]
    fasta // path(fasta)


    main:

    ch_versions            = Channel.empty()
    ch_reports             = Channel.empty()
    ch_empty               = Channel.empty()
    ch_empty_file          = Channel.of("").collect()

    SAMTOOLS_SORT ( bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    PICARD_MARKDUPLICATES( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions.first())
    ch_reports  = ch_reports.mix(PICARD_MARKDUPLICATES.out.metrics)

    SAMTOOLS_INDEX( PICARD_MARKDUPLICATES.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    PICARD_MARKDUPLICATES.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ extract_mapped_input }

    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai)
        .set{ extract_umap_input }

    EXTRACT_MAPPED( extract_mapped_input, fasta )
    ch_versions = ch_versions.mix(EXTRACT_MAPPED.out.versions.first())

    EXTRACT_UMAP( bam )
    ch_versions = ch_versions.mix(EXTRACT_UMAP.out.versions.first())

    SAMTOOLS_FASTQ ( EXTRACT_UMAP.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    // merge R1 and R2, as YAHA expects single end reads
    PEAR( SAMTOOLS_FASTQ.out.fastq )
    ch_versions = ch_versions.mix(PEAR.out.versions.first())

    // create pear_reads channel with [assembled, unassembled_R1, unassembled_R2]
    PEAR.out.assembled
        .join(PEAR.out.unassembled)
        .map{meta, r1,r2,r3 ->
            def meta_tmp = meta
            meta_tmp["single_end"]= true
            [meta_tmp, [r1,r2,r3]]}
        .set{ reads_to_merge }

    // create channel with [ val(meta), [assembled, unassembled_R1, unassembled_R2] ]
    // ch_reads_to_merge.mix(bam.map{meta, bam -> meta})
    // ch_reads_to_merge.mix(pear_reads)

    CAT_FASTQ( reads_to_merge )

    emit:
    aligned         = EXTRACT_MAPPED.out.bam //EXTRACT_MAPPED.out.bam // channel: [ val(meta), path(merged fastq) ]
    unaligned_reads = CAT_FASTQ.out.reads   // channel: [ val(meta), path(merged fastq) ]
    reports         = ch_reports
    versions        = ch_versions           // channel: [ versions.yml ]
}
