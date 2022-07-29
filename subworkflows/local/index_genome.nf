

//
// index the genome in various ways for various tools
//

include { SAMTOOLS_FAIDX } from "${projectDir}/modules/nf-core/modules/samtools/faidx/main"
include { BWAMEM2_INDEX  } from "${projectDir}/modules/nf-core/modules/bwamem2/index/main"
include { YAHA_INDEX     } from "${projectDir}/modules/local/yaha/index/main"

workflow INDEX_GENOME {
    take:
    aligners
    fasta // [val(meta), path(fasta)]

    main:

    ch_versions      = Channel.empty()
    ch_reports       = Channel.empty()

    def fasta_meta = ["id": "genome"]
    ch_fasta_with_meta = Channel.from( fasta_meta )
                                .combine(fasta)

    // index genome with samtools faidx if not input
    if(params.fasta_fai){
        ch_fai = Channel.fromPath(params.fasta_fai).collect()
    } else{
        SAMTOOLS_FAIDX ( ch_fasta_with_meta )
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())
        ch_fai = SAMTOOLS_FAIDX.out.fai.collect()
    }

    //
    // index the genome with bwamem2 index if not passed
    //
    if(params.bwamem2_index && aligners.split(',').contains('bwamem2')){
        ch_bwamem2_index = Channel.fromPath(params.bwamem2_index).collect()
    } else if(aligners && aligners.split(',').contains('bwamem2')){
        BWAMEM2_INDEX ( fasta )
        ch_bwamem2_index = BWAMEM2_INDEX.out.index.collect()
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    } else{
        ch_bwamem2_index = Channel.empty()
    }

    if(params.yaha_index && params.yaha_nib2 && aligners.split(',').contains('ngm_yaha')){
        ch_yaha_index = Channel.fromPath(params.yaha_index).collect()
        ch_yaha_nib2  = Channel.fromPath(params.yaha_nib2).collect()
    } else if(aligners && aligners.split(',').contains('ngm_yaha')){
        YAHA_INDEX( ch_fasta_with_meta )
        ch_yaha_index = YAHA_INDEX.out.index.collect()
        ch_yaha_nib2  = YAHA_INDEX.out.nib2.collect()
        ch_versions   = ch_versions.mix(YAHA_INDEX.out.versions.first())
    } else{
        ch_yaha_index = Channel.empty()
        ch_yaha_nib2  = Channel.empty()
    }

    emit:
    fai           = ch_fai     // channel: [val(meta), path(fai)]
    bwamem2_index = ch_bwamem2_index // path(bwamem2 directory)
    yaha_index    = ch_yaha_index    // channel: [ val(meta), path(yaha index file -- the one with extension .X...)]
    yaha_nib2     = ch_yaha_nib2     // channel: [ val(meta), path(yaha genome .nib2 file) ]
    versions      = ch_versions      // channel: [ versions.yml ]
}


