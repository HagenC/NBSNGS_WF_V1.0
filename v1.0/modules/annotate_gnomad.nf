
process ANNOTATE_GNOMAD {
tag "Annotate gnomad v.4.1 - $sampleId - $flowcell"

//publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.LA.CLN.GNOMAD.vcf", mode:'copy'

input:
 tuple val(sampleId), val(flowcell), path(vcf_file) , path(vcf_index)
 
 output:
 tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.LA.CLN.GNOMAD.vcf"),  emit: CLNGNOMAD_la_vcf_file

script:
"""
	SnpSift annotate ${params.gnomad} $vcf_file      			                     \
	-info  "AF_nfe"  > ${sampleId}.GATK.LA.CLN.GNOMAD.vcf  

	
"""

}
