#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]


#nxf_work_dir <-  "/srv/data/ILMN/PIPELINES/v2_NbsNgsWF"
source(paste0(nxf_work_dir,"/bin/dependencies.R"))

#Importing sampleheets: ----
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/SampleSheets.sqlite"))
sampleProject <- c("NBS-NGS", "nbs-ngs", "NBS_NGS") %>% paste(collapse = "', '") %>%   paste0("'", ., "'") 
query_opl <- paste0(" SELECT * FROM SAMPLESHEETS WHERE Sample_Project IN (",sampleProject,") ")
NBSNGS_sampleSheet <- dbGetQuery(con, query_opl) %>% select(SampleID_Flowcell, Description)
dbDisconnect(con) 

#Importing VOIS: ----
#Importing phenotype-gene list
Phenotype.gene.list <- fread(paste0(nxf_work_dir, "/assets/phenotype/phenotype.gene.list")) %>% pull(Gene.ID)
GENE_filter_SQLformat <- paste0("'", paste(Phenotype.gene.list, collapse = "','"), "'")
#Founder vars
if(file.exists(paste0(nxf_work_dir,"/assets/phenotype/FounderVariants.bed"))== "TRUE"){
  FOUNDER_VARIANTS <- fread(paste0(nxf_work_dir, "/assets/phenotype/FounderVariants.bed"))
}
effect_filter_SQLformat <- paste0("'", paste(keep.SnpEff.EFFECTs, collapse = "','"), "'")

#Extract OPLs
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/variants.sqlite"))
query_opl <- paste0("
   SELECT *  FROM OPL
   WHERE   (EFFECT IN (", effect_filter_SQLformat, ") OR CLNSIG LIKE '%pathogen%')
   AND CLNSIG NOT IN ('Benign', 'Benign/Likely_benign', 'Likely_benign')
   AND GENE in (",GENE_filter_SQLformat,")
 ")

OPL_VOIS <- dbGetQuery(con, query_opl) %>% left_join(., NBSNGS_sampleSheet, by = "SampleID_Flowcell") %>% replace_na(list(AF_nfe = 0))  %>%
  mutate(Pathogenic = grepl( "Pathogenic|Likely_pathogenic", CLNSIG)) %>% filter(AF_nfe < 0.01 | Pathogenic == "TRUE")  %>% select(-Pathogenic)


VOIS_extraction_list <- OPL_VOIS %>% mutate(CHR_POS = paste0(CHROM, "-", POS)) %>% pull(CHR_POS)
VOIS_extraction_list_SQLformat <- paste0("'", paste(VOIS_extraction_list, collapse = "','"), "'")

#Extract COV
query_coverage <- paste0("
  SELECT *
  FROM COV
  WHERE  (CHR_POS IN (", VOIS_extraction_list_SQLformat, "))

")
coverage_VOIS <- dbGetQuery(con,query_coverage)
dbDisconnect(con) 

#Combning extracts
coverage_VOIS_count <- coverage_VOIS %>% group_by(CHR_POS) %>% summarise(PASSED.alleles = length(which( DP >= 10))*2, FAILED.alleles = length(which( DP < 10))*2)
VOI_DB_GENE <- OPL_VOIS %>%   unite(VARID, c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>% select(VARID, GENE) %>% filter(!duplicated(VARID)) 
if(nrow(OPL_VOIS) > 0){
VOI_DB <- OPL_VOIS %>%   unite(VARID, c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>% group_by(VARID) %>%
  mutate(GT = ifelse(GT == "0|1", "0/1", GT)) %>% mutate(GT = ifelse(GT == "1|1", "1/1", GT)) %>%  count(GT) %>% 
  spread(GT, n) 
VOI_DB <- VOI_DB %>%
  mutate(`0/1` = ifelse(!("0/1" %in% names(VOI_DB)), 0, ifelse(is.na(`0/1`), 0, `0/1`)),
         `1/1` = ifelse(!("1/1" %in% names(VOI_DB)), 0, ifelse(is.na(`1/1`), 0, `1/1`))) %>%
  mutate(CHR_POS = str_extract(VARID, "chr[[:digit:]]+-[[:digit:]]+|chr[[:alpha:]]-[[:digit:]]+")) %>% replace_na(list(`0/1` = 0, `1/1` = 0)) %>% mutate(alleles = `0/1` + `1/1`*2) %>%  
  left_join(., coverage_VOIS_count, by = "CHR_POS")  %>% mutate(NBSNGS_AF = alleles/PASSED.alleles) %>% select(-alleles, -CHR_POS) %>% left_join(., VOI_DB_GENE, by = "VARID") %>% 
  select(VARID, GENE, everything()) %>% rename( Heterozygote = `0/1`, Homozygote = `1/1`) %>% mutate(NBSNGS_AF = round(NBSNGS_AF, digits = 5))
}else{
  VOI_DB <- data.frame(VARID = NA,    GENE = NA,    Heterozygote = NA,   Homozygote = NA,      PASSED.alleles = NA,  FAILED.alleles = NA,  NBSNGS_AF = NA)
}
# rename(Unknown.multiallelic = `./.`, Heterozygote = `0/1`, Homozygote = `1/1`)


covarage_PHENO <- coverage_VOIS %>% left_join(., NBSNGS_sampleSheet, by = "SampleID_Flowcell") %>% group_by(CHR_POS, Description) %>% 
  summarise(PASSED.alleles = length(which( DP >= 10))*2, FAILED.alleles = length(which( DP < 10))*2, .groups = "keep") %>% unite(CHR_POS_PHENOTYPE,c(CHR_POS, Description), sep = "_" )
if(nrow(OPL_VOIS) > 0){
VOI_DB_PHENOTYPE <- OPL_VOIS %>%   unite(VARID, c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>% group_by(VARID, Description) %>% 
  mutate(GT = ifelse(GT == "0|1", "0/1", GT)) %>% mutate(GT = ifelse(GT == "1|1", "1/1", GT)) %>%  count(GT) %>% 
  spread(GT, n) 
VOI_DB_PHENOTYPE <- VOI_DB_PHENOTYPE %>%
  mutate(`0/1` = ifelse(!("0/1" %in% names(VOI_DB_PHENOTYPE)), 0, ifelse(is.na(`0/1`), 0, `0/1`)),
         `1/1` = ifelse(!("1/1" %in% names(VOI_DB_PHENOTYPE)), 0, ifelse(is.na(`1/1`), 0, `1/1`))) %>%
  mutate(CHR_POS = str_extract(VARID, "chr[[:digit:]]+-[[:digit:]]+|chr[[:alpha:]]-[[:digit:]]+")) %>% replace_na(list(`0/1` = 0, `1/1` = 0))  %>%
  mutate(alleles = `0/1` + `1/1`*2) %>% 
  unite(CHR_POS_PHENOTYPE,c(CHR_POS, Description), sep = "#" ) %>% 
  left_join(., covarage_PHENO, by = "CHR_POS_PHENOTYPE") %>% separate(CHR_POS_PHENOTYPE, c("CHR_POS", "Phenotype"), sep = "#") %>% mutate(Phenotype_NBSNGS_AF = alleles/PASSED.alleles) %>% select(-alleles, -CHR_POS) %>% left_join(., VOI_DB_GENE, by = "VARID") %>% 
  select(VARID, GENE, everything())  %>% rename(Heterozygote = `0/1`, Homozygote = `1/1`) %>% mutate(Phenotype_NBSNGS_AF = round(Phenotype_NBSNGS_AF, digits = 5))
}else{
  VOI_DB_PHENOTYPE <- data.frame(VARID = NA,    GENE = NA,    Heterozygote = NA,   Homozygote = NA,      PASSED.alleles = NA,  FAILED.alleles = NA,  Phenotype_NBSNGS_AF = NA)
}


#Export 
if(exists("VOI_DB") == "TRUE"){
  fwrite(VOI_DB, "VOI_DB", sep = "\t")
  fwrite(VOI_DB_PHENOTYPE , "VOI_DB_PHENOTYPE", sep = "\t" )
}else{
  VOI_DB <-  data.frame(VARID = NA,  GENE = NA,  Heterozygote = NA , Homozygote = NA,  PASSED.alleles = NA, FAILED.alleles = NA, NBSNGS_AF = NA)
  fwrite(VOI_DB, "VOI_DB", sep = "\t")
  VOI_DB_PHENOTYPE <- data.frame(VARID = NA,  GENE = NA, Phenotype = NA,  Heterozygote = NA , Homozygote = NA,  PASSED.alleles = NA, FAILED.alleles = NA, Phenotype_NBSNGS_AF = NA)
  fwrite(VOI_DB_PHENOTYPE , "VOI_DB_PHENOTYPE", sep = "\t" )
}



