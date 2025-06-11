#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
VOI_DB_path <- args[2]
VOI_DB_PHENOTYPE_path <- args[3]
version = args[4]

#nxf_work_dir <-  "/srv/data/ILMN/PIPELINES/NBSNGS_WF_v1.0/v1.0/"
#version <- "v.1.0"
#VOI_DB_path <- "/srv/data/ILMN/PIPELINES/NBSNGS_WF_v1.0/v1.0/assets/temp/VOI_DB"
source(paste0(nxf_work_dir,"/bin/dependencies.R"))



#Importing samplesheet table
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/SampleSheets.sqlite"))
sampleProject <- c("NBS-NGS", "nbs-ngs", "NBS_NGS") %>% paste(collapse = "', '") %>%   paste0("'", ., "'") 
query_samplesheet <- paste0(" SELECT * FROM SAMPLESHEETS WHERE Sample_Project IN (",sampleProject,")")
samplesheet_NBSNGS <- dbGetQuery(con, query_samplesheet) %>% select(SampleID_Flowcell, Description, Flowcell)
dbDisconnect(con) 

#Importing OPL SampleID_Flowcell table
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/variants.sqlite"))
OPL_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM OPL")
dbDisconnect(con) 

DONE_REPORTS <- data.frame(REPORT = list.files(paste0(nxf_work_dir, "/../REPORTS_",version,"/"), "*report.html$", recursive = TRUE)) %>% 
  separate(REPORT, c("Flowcell", "html"), sep ="/") %>% separate(html , c("SampleID", "Flowcell.sanity", "Phenotype"), sep ="-") %>% 
  unite(SampleID_Flowcell, c(SampleID, Flowcell), sep = "-") %>% mutate(Phenotype = gsub(".report.html", "", Phenotype))

SAMPLES_FOR_REPORTING <- samplesheet_NBSNGS %>% filter(SampleID_Flowcell %in% OPL_SampleID_Flowcells$SampleID_Flowcell) %>% filter(!SampleID_Flowcell %in% DONE_REPORTS$SampleID_Flowcell)


#Importing ---- #############################################################################################################################################
#Importing and prepping Clinvar pathogenic subset.  
CLN <- suppressWarnings(fread(paste0(nxf_work_dir, "/assets/temp/clinvarPathogenicTargetsubet.tsv")) %>% separate(GENEINFO, c("Symbol", "Gene.No"), sep = ":") %>% 
 mutate(CHROM = paste0("chr", CHROM)) %>% unite(CHR_POS, CHROM:POS, sep = "-"))
#Importing beds
bedfile_targets <- fread(paste0(nxf_work_dir,"/assets/temp/bedfile_target.bed"))
bedfile_coverage_expanded <- fread(paste0(nxf_work_dir,"/assets/temp/bedfile_coverage_expanded.bed"))
bedfile_call_expanded <- fread(paste0(nxf_work_dir,"/assets/temp/bedfile_call_expanded.bed"))
#Importing phenotype-gene list
Phenotype.gene.list <- fread(paste0(nxf_work_dir, "/assets/phenotype/phenotype.gene.list"))
#Founder vars
if(file.exists(paste0(nxf_work_dir,"/assets/phenotype/FounderVariants.bed"))== "TRUE"){
FOUNDER_VARIANTS <- fread(paste0(nxf_work_dir, "/assets/phenotype/FounderVariants.bed"))
}
VOI_DB_import <- fread(VOI_DB_path) %>% select(-GENE)

#i <- 33
print(paste0("Processing ", nrow(SAMPLES_FOR_REPORTING), " samples:"))
if(nrow(SAMPLES_FOR_REPORTING) > 0){
for(i in 1:nrow(SAMPLES_FOR_REPORTING)){
targetSampleID <- SAMPLES_FOR_REPORTING$SampleID_Flowcell[i]
targetPhenotype <-   SAMPLES_FOR_REPORTING$Description[i]
targetFlowcell <-   SAMPLES_FOR_REPORTING$Flowcell[i]
#target genes: ---- #############################################################################################################################################
PHENOTYPE.GENES <- Phenotype.gene.list %>% filter(Phenotype == targetPhenotype) %>% pull(Gene.ID)
if(length(PHENOTYPE.GENES) == 0){
  PHENOTYPE.GENES <- Phenotype.gene.list$Gene.ID
}


#Importing SampleID_Flowcell data : ---- #####################################################################################################################
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/variants.sqlite"))
#VOIS table ----
#effect_filter_SQLformat <- paste0("'", paste(keep.SnpEff.EFFECTs, collapse = "','"), "'")
effect_like_sql<- paste0("EFFECT LIKE '%", keep.SnpEff.EFFECTs, "%'", collapse = " OR ")
GENE_filter_SQLformat <- paste0("'", paste(PHENOTYPE.GENES, collapse = "','"), "'")

 query_opl <- paste0("SELECT *  FROM OPL WHERE SampleID_Flowcell =  '",targetSampleID,"'")
 query_opl <- paste0("
    SELECT *  FROM OPL WHERE SampleID_Flowcell =  '",targetSampleID,"'
    AND ((", effect_like_sql, ") OR CLNSIG LIKE '%pathogen%')
    AND (CLNSIG NOT IN ('Benign', 'Benign/Likely_benign', 'Likely_benign') OR CLNSIG IS NULL)
    AND GENE in (",GENE_filter_SQLformat,") ")
OPL_VOIS <- dbGetQuery(con, query_opl)  %>% unite(VARID, c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>% replace_na(list(AF_nfe = 0))  %>%
  mutate(Pathogenic = grepl( "Pathogenic|Likely_pathogenic", CLNSIG)) %>% filter(AF_nfe < 0.01 | Pathogenic == "TRUE")   %>% select(-Pathogenic)


  
#RAW - sampledID
query_RAW <- paste0(" SELECT * FROM RAW WHERE SampleID_Flowcell = '",targetSampleID,"'")
RAW <- dbGetQuery(con,query_RAW) %>% select(VARID, SB, DP4) %>% filter(!duplicated(VARID))
query_LOFREQ <- paste0(" SELECT * FROM LOFREQ WHERE SampleID_Flowcell = '",targetSampleID,"'")
LOFREQ  <- dbGetQuery(con,query_LOFREQ)

if(nrow(OPL_VOIS) != 0){
VOIS_targetSampleID_annotated <- OPL_VOIS %>% left_join(., RAW, by = "VARID") %>% mutate(LOFREQ.supported = ifelse(VARID %in% LOFREQ$VARID, "Supported", "not-supported")) %>%
  left_join(., VOI_DB_import, by = "VARID") %>%
  select(VARID, RS, AF, SB, DP4, GENE, EFFECT, HGVS_C, HGVS_P, CLNSIG, AF_nfe, NBSNGS_AF, LOFREQ.supported, SampleID_Flowcell) %>%
  mutate(AF = round(AF, digits = 5), NBSNGS_AF = round(NBSNGS_AF, digits = 5))
}else{
  VOIS_targetSampleID_annotated  <- data.frame(matrix(rep("-", times = length(names(OPL_VOIS))),
                                      nrow =1,  dimnames = list(1, names(OPL_VOIS))))   %>% mutate(SB = "-", DP4 = "-", NBSNGS_AF = "-", LOFREQ.supported = "-") %>%
    select(VARID, RS, AF, SB, DP4, GENE, EFFECT, HGVS_C, HGVS_P, CLNSIG, AF_nfe, NBSNGS_AF, LOFREQ.supported, SampleID_Flowcell)
}

#Target region variation: ----
query_variation <- paste0("
   SELECT *  FROM OPL WHERE SampleID_Flowcell =  '",targetSampleID,"'
   AND GENE in (",GENE_filter_SQLformat,") ")
VARIATION_targetSampleID <- dbGetQuery(con, query_variation)  %>% unite(VARID, c(CHROM, POS, REF, ALT), sep = "-", remove = FALSE) %>%  replace_na(list(AF_nfe = 0)) %>% 
  select(VARID, RSID, AF,AF_nfe, GENE, EFFECT, HGVS_C, HGVS_P, CLNSIG) %>% mutate(AF = round(AF, digits = 5))
if(nrow(VARIATION_targetSampleID )==0){
  VARIATION_targetSampleID  <- data.frame(matrix(rep("-", times = length(names(VARIATION_targetSampleID))),
                                                      nrow =1,  dimnames = list(1, names(VARIATION_targetSampleID))))
}


#VOI Phenotype count

VOI_DB_PHENOTYPE_import <- fread(VOI_DB_PHENOTYPE_path) %>% filter(VARID %in%VOIS_targetSampleID_annotated$VARID)
if(nrow(VOI_DB_PHENOTYPE_import )==0){
  VOI_DB_PHENOTYPE_import  <- data.frame(matrix(rep("-", times = length(names(VOI_DB_PHENOTYPE_import))),
                                                 nrow =1,  dimnames = list(1, names(VOI_DB_PHENOTYPE_import))))
}

#Coverage - sampledID ----
query_coverage <- paste0("
  SELECT * FROM COV WHERE SampleID_Flowcell = '",targetSampleID,"'")
coverage <- dbGetQuery(con,query_coverage)
#Coverage sampleID exon
Coverage_targetSampleID_genesExons  <- coverage%>% inner_join(., bedfile_targets , by = c("CHR_POS" = "chr.pos")) %>%
  mutate(Chromosome = str_extract(CHR_POS, "chr[[:digit:]]+|chrX|chrY")) %>%
  unite(Region, c(Chromosome,Start,End), sep = "-") %>% filter(GeneID %in% PHENOTYPE.GENES) %>% 
  group_by(GeneID, Exon, Region) %>%
  summarise(meanCoverage = round(mean(DP), digits = 1), size = n(), covered = length(which(DP >= 10)), .groups = "keep") %>% mutate(coveredPercent = round(covered/size*100, digits = 1)) %>%
  select(-size, -covered) %>% rename(GENE = GeneID)
#SEX:
query_sex <- paste0("
   SELECT *  FROM SEX WHERE SampleID_Flowcell =  '",targetSampleID,"'" )
SEX_targetSampleID <- dbGetQuery(con, query_sex) %>% pull(Predicted_Sex)
if(length(SEX_targetSampleID) == 0){
  SEX_targetSampleID  <- "Unknown"
}

query_mtDNA <- paste0("
   SELECT *  FROM mtDNA WHERE SampleID_Flowcell =  '",targetSampleID,"'" )
mtDNA_haplogroup <- dbGetQuery(con, query_mtDNA)  %>% mutate(Qual = paste0(Haplogroup," (Quality=",Quality,";n=", Variants, ")")) %>% 
  select(-id, -SampleID_Flowcell, -Quality, -Variants, -Haplogroup) %>% .[1,]
if(length(mtDNA_haplogroup) == 0){
  mtDNA_haplogroup  <- "Insufficient data"
}

query_ANCESTRY <- paste0("
   SELECT *  FROM ANCESTRY WHERE SampleID_Flowcell =  '",targetSampleID,"'" )
Ancestry_predict <- dbGetQuery(con, query_ANCESTRY)  %>%  select(INFO) 
if(nrow(Ancestry_predict) == 0){
  Ancestry_predict  <- "Insufficient data (Probability=0;Info=0;Accuracy=0)"
}




dbDisconnect(con)

#CLINVAR
if(length(PHENOTYPE.GENES) > 20){
clinvar.tagets <- CLN %>% filter(Symbol %in% PHENOTYPE.GENES)
cov.clinvar   <- left_join(clinvar.tagets,coverage , by = "CHR_POS") %>%  select(-ID) %>% select(CHR_POS,REF, ALT,  DP ,RS, CLNSIG,CLNREVSTAT, Symbol) %>% arrange(DP) %>%
  filter(CHR_POS %in% bedfile_coverage_expanded$chr.pos) %>% unite(Var.ID, CHR_POS:ALT, sep = "-") %>% left_join(., RAW, by = c("Var.ID" = "VARID")) %>%
  replace_na(list(AF.raw = "-", DP4.raw = "-", nt.raw = "-")) %>% filter(DP < 10 | is.na(DP) | !is.na(DP4)) %>% rename(VARID = Var.ID, GENE = Symbol)
}else{
  clinvar.tagets <- CLN %>% filter(Symbol %in% PHENOTYPE.GENES)
  cov.clinvar   <- left_join(clinvar.tagets,coverage , by = "CHR_POS") %>%  select(-ID) %>% select(CHR_POS,REF, ALT,  DP ,RS, CLNSIG,CLNREVSTAT, Symbol) %>% arrange(DP) %>%
    filter(CHR_POS %in% bedfile_coverage_expanded$chr.pos) %>% unite(Var.ID, CHR_POS:ALT, sep = "-") %>% left_join(., RAW, by = c("Var.ID" = "VARID")) %>%
    replace_na(list(AF.raw = "-", DP4.raw = "-", nt.raw = "-"))%>% rename(VARID = Var.ID, GENE = Symbol) 
  
}
#FOUNDER
n.founders <- FOUNDER_VARIANTS %>% filter(GENE %in% PHENOTYPE.GENES)
if(nrow(n.founders) > 0){
FOUNDER_TABLE <- FOUNDER_VARIANTS %>% left_join(.,VARIATION_targetSampleID, by = "VARID" ) %>% left_join(., RAW, by = "VARID") %>%
  left_join(., coverage, by = "CHR_POS") %>% select(INFO, GENE.x, VARID, AF, DP4, DP ) %>% rename(GENE = GENE.x)
}else{
  FOUNDER_TABLE  <- data.frame(INFO = "-", GENE = "-", VARID = "-", AF = "-", DP4 = "-", DP = "-")
}



#Rendering ----
phenotype_suffix <- targetPhenotype %>% gsub("/","_", .)
output_file <- paste0(targetSampleID,"-", phenotype_suffix,".report.html")

# Define the parameters for the report
params_list = list(
  sampleID = targetSampleID,
  flowcellID = targetFlowcell,
  phenotype= targetPhenotype,
  screenedGenes= PHENOTYPE.GENES ,
  targetVoisAF= VOIS_targetSampleID_annotated,
  targetVariation= VARIATION_targetSampleID,
  targetVoisPopulationAF= VOI_DB_PHENOTYPE_import,
  geneExonCov= Coverage_targetSampleID_genesExons,
  CPR= "-",
  covClinvarPathogenic= cov.clinvar ,
  founderCovTargets= FOUNDER_TABLE, 
  Sex = SEX_targetSampleID, 
  mtDNA = mtDNA_haplogroup,
  ancestry = Ancestry_predict,
  Version = version
)

# Render the report
rmarkdown::render(input =paste0(nxf_work_dir,"/bin/report.Rmd"),
       output_file = output_file,
       output_dir = paste0(nxf_work_dir,"/../REPORTS_",version,"/",targetFlowcell),
       params = params_list)

}
}
fwrite(SAMPLES_FOR_REPORTING , "Processed_SampleIDs.txt", sep = "\t")

# ##fwrite(cov.clinvar, paste0("./NBS.NGS.output/TMP/", ID ,".cov.clinvar" ), sep = "\t") 
# #EXPORT: ---- #############################################################################################################################################
# fwrite(TARGET_VOIS.AF, paste0(SampleId,".targetVoisAF"), sep = "\t")
# fwrite(cov.clinvar  , paste0(SampleId, ".covClinvarPathogenic"), sep = "\t")
# fwrite(founder.targets.cov , paste0(SampleId,".founderCovTargets"), sep = "\t")
# fwrite(GeneExonCov, paste0(SampleId,".geneExonCov"), sep = "\t")
# 
# fwrite(TARGET_VARIATION.out, paste0(SampleId,".targetVariation"), sep = "\t")
# fwrite(data.frame(Genes = PHENOTYPE.GENES), paste0(SampleId,".screenedGenes"))
# fwrite(voiddbPhenotypeI, paste0(SampleId,".voiDBPhenotypeTargets"), sep = "\t")


