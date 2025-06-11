process APPLY_BQSR {
    tag "Apply score recalibration - $sampleId - $flowcell"

 
    input:
    tuple val(sampleId),val(flowcell), path(bam_file, stageAs: "raw/*"), path(bqsr_table)
    

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.bam") , emit: corrected_bam_file

    script:
    """
    gatk ApplyBQSR        \
        -R ${params.reference_genome}    \
        -I ${bam_file}                   \
        -bqsr ${bqsr_table}              \
        -O "${sampleId}.bam" 
    """
}
