process BQSR {
    tag "BQSR - $sampleId - $flowcell"

    input:
    tuple val(sampleId), val(flowcell), path(bam_file)
   

    output:
    tuple val(sampleId),val(flowcell) ,path(bam_file), path("${sampleId}.BQSR.table"), emit: bqsr_file

    script:
    """
    gatk BaseRecalibrator \
        -I ${bam_file} \
        -R ${params.reference_genome} \
        --known-sites ${params.mills} \
        --known-sites ${params.g1000} \
        -O ${sampleId}.BQSR.table
    """

}
