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
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
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
include { PREPARE_INTERVALS        } from "${projectDir}/subworkflows/local/prepare_intervals"
include { ALIGN                    } from "${projectDir}/subworkflows/local/align"
include { CALL_INDIVIDUAL_VARIANTS } from "${projectDir}/subworkflows/local/call_individual_variants"
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
fasta_fai = params.fasta_fai ?
            Channel.fromPath(params.fasta_fai).collect() :
            Channel.empty()

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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_reports  = ch_reports.mix(FASTQC.out.zip.collect{meta, logs -> logs})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // Build intervals if needed
    PREPARE_INTERVALS(
        fasta_fai
    )
    ch_versions = ch_versions.mix(PREPARE_INTERVALS.out.versions.first())
    ch_intervals            = PREPARE_INTERVALS.out.intervals_bed
    ch_intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi
    ch_intervals_combined   = PREPARE_INTERVALS.out.intervals_bed_combined

    ALIGN (
        INPUT_CHECK.out.reads,
        fasta,
        fasta_fai,
        ch_intervals_combined
    )
    ch_reports  = ch_reports.mix(ALIGN.out.report)
    ch_versions = ch_versions.mix(ALIGN.out.versions.first())

    // CALL_INDIVIDUAL_VARIANTS (
    //     ALIGN.out.bam_bai,
    //     fasta,
    //     fasta_fai,
    //     [], // bwa index for tiddit not used
    //     ch_intervals,
    //     ch_intervals_bed_gz_tbi,
    //     ch_intervals_combined
    // )
    // ch_versions = ch_versions.mix(CALL_INDIVIDUAL_VARIANTS.out.versions.first())

    // ANNOTATE (
    //     CALL_INDIVIDUAL_VARIANTS.out.freebayes_vcf
    // )
    // ch_reports  = ch_reports.mix(ANNOATE.out.reports)
    // ch_versions = ch_versions.mix(ANNOTATE.out.versions.first())

    // if(params.bsa2){
    //     BSA2 (
    //         ANNOTATE.out.freebayes.collect()
    //     )
    //     ch_versions = ch_versions.mix(BSA2.out.versions.first())
    // }

    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

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
