process HAPLOTYPECALLER {
    tag "HaplotypeCaller - $sampleId - $flowcell"
  
    publishDir "${params.outputDir}/../log_${params.version}/${flowcell}/haplotypecaller/", pattern: "${sampleId}.haplotypecaller.log", mode:'copy'
    publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.g.vcf.gz", mode:'copy'
    publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.g.vcf.gz.tbi", mode: 'copy'

    input:
    tuple val(sampleId), val(flowcell) , path(bam_file), path(bam_index_file)
  

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.g.vcf.gz"), path("${sampleId}.GATK.g.vcf.gz.tbi"), emit: vcf_file
    path "${sampleId}.haplotypecaller.log"

    script:
    """
    gatk --java-options "-Xmx120g -XX:-UsePerfData"       \
        HaplotypeCaller                                 \
        -R ${params.reference_genome}                   \
        -I ${bam_file}                                  \
        -L ${params.call_bed}                       \
        -O "${sampleId}.GATK.g.vcf.gz"                   \
        --pair-hmm-implementation FASTEST_AVAILABLE \
        --native-pair-hmm-threads ${task.cpus}          \
        --max-alternate-alleles 3                       \
        -ERC GVCF                                       \
        > >(tee "${sampleId}.haplotypecaller.log") 2>&1 
    """
  
}
