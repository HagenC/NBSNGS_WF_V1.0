process INSERT_SIZE_METRICS {
    tag "$sampleId"
    // Collect insert size metrics for bam file


    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/insert_size_metrics/", pattern: "${sampleId}.insert.txt", mode:'copy'
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/insert_size_metrics/", pattern: "${sampleId}.pdf", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file)
    

    output:
    path("${sampleId}.insert.txt"), emit: metrics_file

    script:
    
    """
    gatk CollectInsertSizeMetrics           \
        -I ${bam_file}                      \
        -O ${sampleId}.insert.txt                  \
        -H ${sampleId}.pdf
    """

}
