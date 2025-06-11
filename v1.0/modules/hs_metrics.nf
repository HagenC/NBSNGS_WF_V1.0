process HSMETRICS {
    tag "$sampleId - $flowcell"
   
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/hs_metrics/", pattern: "${sampleId}.hs.txt", mode:'copy'

    input:
    tuple val(sampleId),val(flowcell), path(bam_file) 
    

    output:
    path("${sampleId}.hs.txt"), emit: metrics_file
    

    script:
    
    """
    gatk  BedToIntervalList \
    I=${params.coverage_bed} \
    O="target_region.dict" \
    SD=${params.reference_genome}
           
    gatk CollectHsMetrics    \
        -I ${bam_file}                      \
        -O ${sampleId}.hs.txt                  \
        -R ${params.reference_genome}                        \
        --BAIT_INTERVALS "target_region.dict"      \
        --TARGET_INTERVALS "target_region.dict"   \
	--PER_TARGET_COVERAGE "target_region.dict"
 
    """
  }  
