process CNVPYTOR_VIEW {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cnvpytor=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.2.1--pyhdfd78af_0':
        'quay.io/biocontainers/cnvpytor:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pytor)
    val bin_size
    val output_format
    path config_file
    path gc_content

    output:
    tuple val(meta), path("*.vcf"), emit: vcf      , optional: true
    tuple val(meta), path("*.tsv"), emit: tsv      , optional: true
    tuple val(meta), path("*.xls"), emit: xls      , optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // TODO i coudln't get the output suffix to be a string, so i hard coded 
    // it below. it kept coming out [tsv]
    def output_suffix = "${output_format}" ?: 'tsv'
    def bins   = bin_size ?: '1000'
    def input  = pytor.join(" ")
    def prefix = task.ext.prefix ?: "${meta.id}"
    def conf_arg = config_file ? "-conf ${config_file}" : ''
    """
    cnvpytor \\
    ${conf_arg} \\
    -root ${pytor} \\
    -view ${bin_size} \\
    <<ENDL 
    set print_filename ${prefix}.tsv
    print calls
    ENDL

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/CNVpytor //' )
    END_VERSIONS
    """

    stub:
    def output_suffix = "${output_format}" ?: 'vcf'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${output_suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/CNVpytor //' )
    END_VERSIONS
    """
}