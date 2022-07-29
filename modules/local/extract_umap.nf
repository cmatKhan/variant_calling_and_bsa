//
// Extract unmapped read pairs.
// See http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/
//

process EXTRACT_UMAP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::samtools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.15.1--h1170115_0' :
        'quay.io/biocontainers/samtools:1.15.1--h1170115_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_umap.bam") , emit: bam , optional: true
    tuple val(meta), path("*_umap.cram"), emit: cram, optional: true
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"

    if(meta.single_end){
        """
        samtools view \\
            --threads ${task.cpus-1} \\
            -uf 4 \\
            $args \\
            $input > ${prefix}.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    } else{
        """
        # extract unmapped read whose mate is mapped
        samtools view \\
            --threads ${task.cpus-1} \\
            -u  \\
            -f 4 \\
            -F 264 \\
            ${input}  > ${prefix}_tmp1.bam

        # extract mapped read whose mate is unmapped
        samtools view \\
            --threads ${task.cpus-1} \\
            -u \\
            -f 8 \\
            -F 260 \\
            ${input}  > ${prefix}_tmp2.bam

        # extract reads where both r1 and r2 are unmapped
        samtools view \\
            --threads ${task.cpus-1} \\
            -u \\
            -f 12 \\
            -F 256 \\
            ${input} > ${prefix}_tmp3.bam

        samtools merge \\
            -u \\
            - \\
            ${prefix}_tmp[123].bam | \\
        samtools sort \\
            $args2 \\
            -@ ${task.cpus-1} \\
            -n \\
            - \\
            -o ${prefix}.bam \\
            -T $prefix

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

