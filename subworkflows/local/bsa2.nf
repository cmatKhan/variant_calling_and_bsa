
workflow BSA2 {
    take:
    reads            // channel: [ val(meta), [ reads ] ]
    genome           // channel: file(ref_genome)

    main:

    ch_versions = Channel.empty()
    ch_bam      = Channel.empty()

    if(params.aligner == 'bwamem2') {
        BWAMEM2_ALIGNER (
            reads,
            genome
        )
        ch_bam      = ch_bam.mix(BWAMEM2_ALIGNER.out.bam)
        ch_versions = ch_versions.mix(BWAMEM2_ALIGNER.out.versions)
    } else {
        exit 1, "No aligner specified in params OR aligner: ${params.aligner} is not recognized. "
    }


    emit:
    bam       = ch_bam           // channel: [ val(meta), [ bam ] ]
    versions  = ch_versions      // channel: [ versions.yml ]
}
