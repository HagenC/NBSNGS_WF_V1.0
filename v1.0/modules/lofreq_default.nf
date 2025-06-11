process LOFREQCALL {
    tag "Lofreq default - $sampleId -  $flowcell"
    //source ${params.activate} LOFREQ 
   

    publishDir "${params.outputDir}/../log_${params.version}/${flowcell}/lofreq/", pattern: "${sampleId}.lofreqDefault.log", mode:'copy'
    publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}", pattern: "${sampleId}.lofreqDefault.vcf.gz", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file), path(bai_file)
    
    output:
    tuple val(sampleId), val(flowcell) ,path("${sampleId}.lofreqDefault.vcf.gz"), emit: vcf_file
    path "${sampleId}.lofreqDefault.log"

script:
    
    """
    lofreq call                            \
        --force-overwrite                           \
        --call-indels                               \
        -f ${params.reference_genome}               \
        -l ${params.call_bed}                 \
        -o "${sampleId}.lofreqDefault.vcf.gz"       \
        ${bam_file}                                 \
        > >(tee "${sampleId}.lofreqDefault.log") 2>&1
    """
  
}
