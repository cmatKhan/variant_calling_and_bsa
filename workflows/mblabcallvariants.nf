/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMblabcallvariants.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta, params.dusted_bed ]
for (param in checkPathParamList) {
    if (param) {
        file(param, checkIfExists: true)
    }
}

// Check mandatory parameters
if (params.input) {
    ch_input = file(params.input)
} else {
    exit 1, 'Input samplesheet not specified!'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml",
                                checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ?
                            Channel.fromPath(params.multiqc_config) :
                            Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK              } from "${projectDir}/subworkflows/local/input_check"
include { PREPARE_GENOME             } from "${projectDir}/subworkflows/local/prepare_genome"
include { PREPARE_INTERVALS        } from "${projectDir}/subworkflows/local/prepare_intervals"
include { ALIGN                    } from "${projectDir}/subworkflows/local/align"
include { CALL_INDIVIDUAL_VARIANTS } from "${projectDir}/subworkflows/local/call_individual_variants"
include { CALL_STRUCTURAL_VARIANTS } from "${projectDir}/subworkflows/local/call_structural_variants"
include { CALL_BATCH_VARIANTS } from "${projectDir}/subworkflows/local/call_batch_variants"
include { ANNOTATE                 } from "${projectDir}/subworkflows/local/annotate_individual_variants"
include { BSA2                     } from "${projectDir}/subworkflows/local/bsa2"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from "${projectDir}/modules/nf-core/modules/fastqc/main"
include { MULTIQC                     } from "${projectDir}/modules/nf-core/modules/multiqc/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "${projectDir}/modules/nf-core/modules/custom/dumpsoftwareversions/main"

fasta     = params.fasta     ?
            Channel.fromPath(params.fasta).collect()     :
            Channel.empty()

cnvpytor_genome_conf_ch = params.cnvpytor_conf_file ?
                          Channel.fromPath(params.cnvpytor_conf_file).collect():
                          Channel.empty()

cnv_histogram_bin_size = Channel.from(params.cnv_histogram_bin_size).collect()
cnv_output_format      = Channel.from(params.cnv_output_format).collect()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MBLABCALLVARIANTS {

    // store software versions from all processes
    ch_versions             = Channel.empty()
    // store reports that will be consumed by multiQC
    ch_reports              = Channel.empty()

    annotate_variants_input = Channel.empty()

    ch_dusted_bed           = params.dusted_bed ?
                              Channel.fromPath(params.dusted_bed).collect() :
                              Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    PREPARE_GENOME (
        params.aligners,
        fasta
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Build intervals if needed
    // note that this comes directory from nf-core/sarek. output is just
    // a bed file of the chr in the fasta. It is input to the variant calling
    // and variant annotation steps
    PREPARE_INTERVALS(
        PREPARE_GENOME.out.fai
    )
    ch_versions             = ch_versions.mix(PREPARE_INTERVALS.out.versions.first())
    ch_intervals            = PREPARE_INTERVALS.out.intervals_bed
    ch_intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
    ch_intervals_combined   = PREPARE_INTERVALS.out.intervals_bed_combined

    // Align the reads
    ALIGN (
        params.aligners,
        INPUT_CHECK.out.reads,
        fasta,
        PREPARE_GENOME.out.fai,
        PREPARE_GENOME.out.bwamem2_index,
        PREPARE_GENOME.out.yaha_index,
        PREPARE_GENOME.out.yaha_nib2
    )
    ch_reports  = ch_reports.mix(ALIGN.out.report)
    ch_versions = ch_versions.mix(ALIGN.out.versions.first())

    // Create VCF files -- note this is on individual alignment files.
    // if there are groups, this is run in addition to the call_batch_variants
    // subworkflow
    if(params.call_individual_variants){
        CALL_INDIVIDUAL_VARIANTS (
            ALIGN.out.bam_bai,
            ch_dusted_bed,
            fasta,
            PREPARE_GENOME.out.fai,
            ch_intervals_combined,
            cnvpytor_genome_conf_ch,
            cnv_histogram_bin_size,
            cnv_output_format
        )
        // combine individual VCF and batch VCF into single channel
        annotate_variants_input = annotate_variants_input
                                    .mix(CALL_INDIVIDUAL_VARIANTS.out.freebayes_vcf)
        ch_versions = ch_versions.mix(CALL_INDIVIDUAL_VARIANTS.out.versions.first())
    }

    CALL_STRUCTURAL_VARIANTS (
        ALIGN.out.bam_bai,
        ch_dusted_bed,
        fasta,
        PREPARE_GENOME.out.bwa_index
    )
    annotate_variants_input = annotate_variants_input.mix(CALL_STRUCTURAL_VARIANTS.out.vcf)
    ch_versions = ch_versions.mix(CALL_STRUCTURAL_VARIANTS.out.versions)

    if(params.call_batch_variants){
        // group the bams by the group key in the metadata map. The grouping is
        // input by the user in the samplesheet.
        // input:  [[id: id1, group: 1, ...], bam, bai],
        //         [[id: id2, group: 2, ...], bam2, bai2],
        //         [[id: id3, group: 2, ...], bam3, bai3]
        //         this is combined with the interval channel (by default, queue of bed files by chromosome)
        // output: [[id: group_1, aligner: bwamem2, interval: chr1_1-end, ...], [bam1], [bai1], chr1_bed]
        //         [[id: group_1, aligner: bwamem2, interval: chr2_1-end, ...], [bam1], [bai1], chr2_bed]
        //         ...
        //         [[id: group_2, aligner: bwamem2, interval: chr1_1-end, ...], [bam2], [bai1], chr1_bed]
        //         [[id: group_2, aligner: bwamem2, interval: chr2_1-end, ...], [bam2], [bai1], chr2_bed]
        //         ...
        ALIGN.out.bam_bai.map{meta, bam, bai ->
            return [meta.group, meta.aligner, bam, bai]}
            .groupTuple(by: [0,1])
            .map{group, aligner, bam_list, bai_list ->
                def meta_tmp = ["id":"group_"+group, "aligner":aligner]
                if(bam_list.size() > 1){
                return [meta_tmp,bam_list,bai_list]
                }
            }
            // extract only the interval bed, not the number of intervals
            .combine(ch_intervals.map{it -> it[0]})
            .map{ meta,bam_list,bai_list,interval ->
                    return [add_interval_to_meta(meta, interval),
                            bam_list, bai_list, interval ]}
            .set{ call_batch_variants_input }

        //Run variant calling on the batched alignment files
        CALL_BATCH_VARIANTS (
            call_batch_variants_input,
            ch_dusted_bed,
            fasta,
            PREPARE_GENOME.out.fai,
            PREPARE_GENOME.out.sequence_dict
        )

        // combine individual VCF and batch VCF into single channel
        annotate_variants_input = annotate_variants_input
                                    .mix(CALL_BATCH_VARIANTS.out.freebayes_vcf)
        ch_versions = ch_versions.mix(CALL_BATCH_VARIANTS.out.versions.first())
    }


    // annotate the VCF files
    ANNOTATE (
        annotate_variants_input,
        fasta
    )
    ch_reports  = ch_reports.mix(ANNOTATE.out.reports)
    ch_versions = ch_versions.mix(ANNOTATE.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // //
    // // MODULE: MultiQC
    // //
    // workflow_summary    = WorkflowMblabcallvariants.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files =  Channel.empty().mix(ch_reports.collect().ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow,
                             params,
                             summary_params,
                             projectDir,
                             log,
                             multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// create a meta map for a item in the split fasta channel
// input: a fasta file path that originates from the splitFasta function
//        for example /path/to/...workdir.../fastafilename.1.fasta.gz
//        note that .1.fasta.gz represents the order in which this record was
//        split -- if the chromosomes are in order, then this is chr1.
// output: a meta map in the format [id: contig_<index>], eg [id: contig_1] for the
//         example above
def add_interval_to_meta(Map meta, bed) {

    def new_meta = [:]

    meta.each{ k,v ->
        new_meta[k] = v}

    new_meta.interval = bed.baseName - "bed"

    return new_meta
}
