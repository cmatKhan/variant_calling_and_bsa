quay.io/biocontainers/erds


process ERDS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::erds=v1.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/erds%3A1.1--pl5.22.0_0' :
        'quay.io/biocontainers/erds' }"

    input:
    tuple val(meta), path(bam), path(vcf)
    path fasta

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def threads = task.cpus < 2 ? 1 : task.cpus -1

        """
        erds_pipeline \\
            -b $bam \\
            -v $vcf \\
            -o . \\
            -r $fasta\\
            -c $threads \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            ERDS: 1.1
        END_VERSIONS
        """

}
