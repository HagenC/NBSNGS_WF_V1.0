process GTOVCF {
    tag "Creating vcf - $sampleId - $flowcell"
  
   // publishDir "${params.outputDir}/../log_${params.version}/${flowcell}/haplotypecaller/", pattern: "${sampleId}.haplotypecaller.log", mode:'copy'
    //publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.vcf.gz", mode:'copy'
    //publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.vcf.gz.tbi", mode: 'copy'

    input:
    tuple val(sampleId), val(flowcell) , path(vcf_file), path(vcf_index_file)
  

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.vcf.gz"), path("${sampleId}.GATK.vcf.gz.tbi"), emit: vcffile
    

    script:
    """
    gatk --java-options "-Xmx20g -XX:-UsePerfData"   GenotypeGVCFs    \
        -R ${params.reference_genome}                   \
        -V ${vcf_file}                                  \
        -O "${sampleId}.GATK.vcf.gz"                   
    
    
    """
  
}
