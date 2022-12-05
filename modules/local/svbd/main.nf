process SVDB {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::svdb==2.5.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/svdb:2.5.0--py39h5371cbf_1' :
        'quay.io/biocontainers/svdb:2.5.0--py37h77a2a36_0' }"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input  = vcfs.join(' ').trim()

    """
    svdb \\
        --merge \\
        --vcf $input \\
        $args > ${prefix}_merged.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svdb: \$(echo \$(svdb --version 2>&1) | sed 's/version:\s*v//g' )
    END_VERSIONS
    """
}
