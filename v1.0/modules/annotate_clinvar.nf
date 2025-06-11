
process ANNOTATE_CLINVAR {
tag "Annotate clinvar - $sampleId - $flowcell"

publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.ANNOTATED.vcf.gz", mode:'copy'

input:
 tuple val(sampleId), val(flowcell), path(vcf_file)
 path clinvar_vcf
 path clinvar_tbi
 
 output:
 tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.ANNOTATED.vcf.gz"),  emit: CLN_la_vcf_file

script:
"""
	SnpSift annotate $clinvar_vcf $vcf_file      			                     \
	-info  "CLNSIG,CLNDN,CLNREVSTAT,GENEINFO,CLNVC,RS"  | bgzip -c > ${sampleId}.GATK.ANNOTATED.vcf.gz   

	
"""

}

