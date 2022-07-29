//
// Index the reference genome and align reads
//

include { BWAMEM2_MEM } from "${projectDir}/modules/nf-core/modules/bwamem2/mem/main"

workflow BWAMEM2_ALIGNER {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    bwamem2_index

    main:

    ch_versions      = Channel.empty()

    //
    // Map reads with bwamem2 mem
    //
    sort_bam = false
    BWAMEM2_MEM ( reads, bwamem2_index, sort_bam )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)

    emit:
    bam            = BWAMEM2_MEM.out.bam          // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions      // channel: [ versions.yml ]
}
