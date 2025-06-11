process OPL {
tag "Wrangling one-per-line - $sampleId - $flowcell"

publishDir "${params.outputDir}/../variants_${params.version}/${flowcell}/", pattern: "${sampleId}.GATK.OPL.vcf.gz", mode:'copy'

input:
 tuple val(sampleId), val(flowcell), path(vcf_file)
 
 output:
 tuple val(sampleId), val(flowcell), path("${sampleId}.GATK.OPL.vcf.gz"),  emit: OPL_vcf_file
 val  (sampleId), emit: wait

script:
"""
genes=\$(paste -sd"|" ${params.geneList})

zcat $vcf_file | ${params.OPL} | SnpSift filter "(ANN[*].GENE =~ '^(\${genes})\$')" |  SnpSift extractFields  - CHROM POS REF ALT ID AF FS DP AF_nfe  "CLNSIG" RS CLNVC CLNDN CLNREVSTAT GENEINFO  "ANN[*].GENEID" "ANN[*].HGVS_P" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].GENE" "ANN[*].HGVS_C" GEN[*] | gzip -c > ${sampleId}.GATK.OPL.vcf.gz  

	
"""

}


