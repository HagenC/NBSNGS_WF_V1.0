process ALIGNMENT_METRICS {
    tag "$sampleId"
    // Collect alignment metrics for bam file


    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/alignment_metrics/", pattern: "${sampleId}.align.txt", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell) , path(bam_file)

    
    output:
    path("${sampleId}.align.txt"), emit: metrics_file

    script:
    """
    gatk CollectAlignmentSummaryMetrics     \
        -I ${bam_file}                      \
        -O ${sampleId}.align.txt                  \
        -R ${params.reference_genome}                            \
    """

}
