process MARKDUPLICATES {
    tag "$sampleId - $flowcell"
    
  
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/duplication_metrics/", pattern: "${sampleId}.txt", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file) 
    

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.marked.bam"), emit: marked_bam_file
    path "${sampleId}.txt", emit: metrics_file

    script:
    
    """
    gatk MarkDuplicates                     \
        --I ${bam_file}                     \
        --M ${sampleId}.txt                 \
        --O ${sampleId}.marked.bam

    
    """


}
