process CNVPYTOR_HISTOGRAM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::cnvpytor=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.2.1--pyhdfd78af_0':
        'quay.io/biocontainers/cnvpytor:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(pytor)
    val bin_sizes
    path config_file
    path gc_content

    output:
    tuple val(meta), path("${pytor.baseName}.pytor")	, emit: pytor
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def bins = bin_sizes ?: '1000'
    def args = task.ext.args ?: ''
    def conf_arg = config_file ? "-conf ${config_file}" : ''
    """
    cnvpytor \\
        ${conf_arg} \\
        -root $pytor \\
        $args \\
        -his $bins

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/CNVpytor //' )
    END_VERSIONS
    """

    stub:
    """
    touch ${pytor.baseName}.pytor

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(echo \$(cnvpytor --version 2>&1) | sed 's/CNVpytor //' )
    END_VERSIONS
    """
}
