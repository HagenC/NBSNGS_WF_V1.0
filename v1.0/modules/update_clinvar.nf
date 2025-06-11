process UPDATE_CLINVAR {
tag "Updating Clinvar (GRCh38)"

publishDir "${params.outputDir}/assets/temp/", pattern: "clinvarPathogenicTargetsubet.tsv", mode:'copy'
//publishDir "${params.outputDir}/assets/temp/", pattern: "clinvar.vcf.gz", mode:'copy'
//publishDir "${params.outputDir}/assets/", pattern: "clinvar.vcf.gz.tbi", mode:'copy'


output:
path "clinvarPathogenicTargetsubet.tsv", emit: clinvar_pathogenic_subset
path "clinvar.vcf.gz",       emit: clinvar_vcf
path "clinvar.vcf.gz.tbi",   emit: clinvar_tbi       

script:
"""
wget -O clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -O clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

#remove "chr" from bed
sed 's/^chr//' ${params.coverage_bed} | cut -f1-3 > chrremoved.bed 

#subseet clinvar
bcftools view -R chrremoved.bed clinvar.vcf.gz -o region_clinvar.vcf

#Make one-per-line
cat region_clinvar.vcf | ${params.OPL} | SnpSift  extractFields  - CHROM POS REF ALT ID  "CLNSIG" RS CLNVC CLNDN CLNREVSTAT GENEINFO  > OPL_clinvar.vcf 

#subset to pathogenic
(head -n 1 OPL_clinvar.vcf && grep -i 'patho' OPL_clinvar.vcf) > clinvarPathogenicTargetsubet.tsv

"""
} 
