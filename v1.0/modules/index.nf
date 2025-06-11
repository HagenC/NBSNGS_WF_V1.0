params.outputDir = "."

process INDEX {
    tag "$sampleId  - $flowcell"
    //publishDir "${params.outputDir}/BAM/${flowcell}/", mode: 'copy'
    
    
    input:
    tuple val(sampleId), val(flowcell), path(bam_file) 

    output:
    tuple val(sampleId), val(flowcell), path(bam_file), path("${bam_file}.bai"), emit: bam_file_w_index

    script:
    """
    samtools index   \
        ${bam_file}
    """
}

