process GATK_LA {
 tag "LeftAligning and trimming $sampleId - $flowcell"
 
 //publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.LA.vcf.gz", mode:'copy'
 //publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.LA.vcf.gz.tbi", mode:'copy'
 
 input:
 tuple val(sampleId), val(flowcell), path(vcf_file), path(vcf_index)
 
 output:
 tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.LA.vcf.gz"), path("${sampleId}.GATK.LA.vcf.gz.tbi"),  emit: la_vcf_file

script:
"""
	gatk LeftAlignAndTrimVariants \
   	-R ${params.reference_genome} 		     \
	-V ${vcf_file} 			             \
	-O ${sampleId}.GATK.LA.vcf.gz    		     \
	-split-multi-allelics
	
"""

}

