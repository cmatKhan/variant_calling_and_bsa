process MERGE_NGM_YAHA {
    tag "$yaha_meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(ngm_meta), path(ngm_align)
    tuple val(yaha_meta), path(yaha_align)
    path fasta
    path fai

    output:
    tuple val(yaha_meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val(yaha_meta), path("${prefix}.cram"), optional:true, emit: cram
    path  "versions.yml"                                  , emit: versions

    when:
    "$ngm_meta.id" == "$yaha_meta.id"

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${yaha_meta.id}"
    def file_type = "bam"
    def reference = fasta ? "--reference ${fasta}" : ""

    """
    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${reference} \\
        ${prefix}.${file_type} \\
        $ngm_align $yaha_align

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${yaha_meta.id}${task.ext.suffix}" : "${yaha_meta.id}"
    def file_type = "bam"
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
