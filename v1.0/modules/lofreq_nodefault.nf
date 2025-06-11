
process LOFREQNODEFAULT {
    tag "Lofreq no-default - $sampleId -  $flowcell"
    
   

    publishDir "${params.outputDir}/../log_${params.version}/${flowcell}/lofreq/", pattern: "${sampleId}.lofreqRaw.log", mode:'copy'
    publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}", pattern: "${sampleId}.lofreqRaw.vcf.gz", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file), path(bai_file)
    
    output:
    tuple val(sampleId), val(flowcell) ,path("${sampleId}.lofreqRaw.vcf.gz"), emit: vcf_file
    path "${sampleId}.lofreqRaw.log"

    script:
    
    """
    lofreq call                            \
        --force-overwrite                           \
         --call-indels                               \
        --no-default-filter                         \
        --sig 1                                     \
        --bonf 1                                    \
        -f ${params.reference_genome}               \
        -l ${params.call_bed}                    \
        -o "${sampleId}.lofreqRaw.vcf.gz"             \
        ${bam_file}                                 \
        > >(tee "${sampleId}.lofreqRaw.log") 2>&1
    """
  
}
