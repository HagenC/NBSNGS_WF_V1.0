process ANNOTATE_ONTOLOGY {
tag "Annotate sequence ontology (SnpEff GRCh38.113 Ensemble - homebrew) - $sampleId - $flowcell"

//publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.LA.CLN.GNOMAD.SO.vcf", mode:'copy'

input:
 tuple val(sampleId), val(flowcell), path(vcf_file)
 
 output:
 tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.LA.CLN.GNOMAD.SO.vcf"),  emit: CLNGNOMADSO_la_vcf_file

script:
"""
    export _JAVA_OPTIONS="-Xmx24G"    
 snpEff -c ${params.snpeff_config} -dataDir ${params.snpeff_data} GRCh38.113  $vcf_file > ${sampleId}.GATK.LA.CLN.GNOMAD.SO.vcf  

	
"""

}



