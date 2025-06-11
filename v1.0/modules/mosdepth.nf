process MOSDEPTH {
    tag "mosdepth - ${sampleId}"

   
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/mosdepth/", pattern: "${sampleId}*.per-base.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/mosdepth/", pattern: "${sampleId}*.regions.bed.gz", mode:'copy'
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/mosdepth/", pattern: "${sampleId}*.mosdepth.global.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/mosdepth/", pattern: "${sampleId}*.mosdepth.region.dist.txt", mode:'copy'
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/mosdepth/", pattern: "${sampleId}*.mosdepth.summary.txt", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell) , path(bam), path(bai)
   

    output:
    path "${sampleId}*.per-base.bed.gz", emit: per_base_depth
    path "${sampleId}*.regions.bed.gz", emit: region_depth
    path "${sampleId}*.mosdepth.global.dist.txt", emit: global_dist
    path "${sampleId}*.mosdepth.region.dist.txt", emit: region_dist
    path "${sampleId}*.mosdepth.summary.txt", emit: summary

    script:
    
    """
    mosdepth                                \
        --threads ${task.cpus}              \
        --by ${params.coverage_bed}                 \
        ${sampleId}                     \
        ${bam}
    """
    
}
