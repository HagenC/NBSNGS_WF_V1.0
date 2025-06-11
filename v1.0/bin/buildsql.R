#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
#bedfile_coverage <- args[2]
#bedfile_call <- args[3]
version <- args[2]

source(paste0(nxf_work_dir,"/bin/dependencies.R"))

#nxf_work_dir <-  "~/TEST_TMP/NBSNGS_test/NBSNGS_WF/v1.0"
#version <- "v.1.0"

#Getting sample sheets:
#open connection:

con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/SampleSheets.sqlite"))
# 
 sampleProject <- c("NBS-NGS", "nbs-ngs", "NBS_NGS") %>% paste(collapse = "', '") %>%   paste0("'", ., "'") 
 query_opl <- paste0("
   SELECT *
   FROM SAMPLESHEETS
   WHERE Sample_Project IN (",sampleProject,")
 ")
 
NBSNGS_sampleSheet <- dbGetQuery(con, query_opl) 
#kill connection  
dbDisconnect(con) 


#Importing files 
bedfile_coverage_expanded <- fread(paste0(nxf_work_dir,"/assets/temp/bedfile_coverage_expanded.bed"))
bedfile_call_expanded <- fread(paste0(nxf_work_dir,"/assets/temp/bedfile_call_expanded.bed"))
phenotypeGenes <- fread(paste0(nxf_work_dir,"/assets/phenotype/phenotype.gene.list"))
#Subsetting bed to relevant screening targets
bedfile_coverage_expanded_targets <- bedfile_coverage_expanded   %>% filter(GeneID %in% phenotypeGenes$Gene.ID )


#FIND files ----
OPLs <- data.frame(OPL = list.files(paste0(nxf_work_dir, "/../variants_",version,"/"), "*.GATK.OPL.vcf.gz", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(OPL,stri_locate_last_fixed(OPL, "/")[,1]+1, nchar(OPL))) %>%
  mutate(SampleID_Flowcell  = sub(".GATK.OPL.vcf.gz", "", SampleID_Flowcell )) %>% filter(SampleID_Flowcell  %in% NBSNGS_sampleSheet$SampleID_Flowcell )

LOFREQs <- data.frame(LOFREQ = list.files(paste0(nxf_work_dir, "/../variants_",version,"/"), "*.lofreqDefault.vcf.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(LOFREQ,stri_locate_last_fixed(LOFREQ, "/")[,1]+1, nchar(LOFREQ))) %>%
  mutate(SampleID_Flowcell = sub(".lofreqDefault.vcf.gz", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)

COVs <- data.frame(COV = list.files(paste0(nxf_work_dir, "/../QC_",version,"/"), "*.cov.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(COV,stri_locate_last_fixed(COV, "/")[,1]+1, nchar(COV))) %>%
  mutate(SampleID_Flowcell = sub(".cov.gz", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)

RAWs <- data.frame(RAW = list.files(paste0(nxf_work_dir, "/../variants_",version,"/"), "*.lofreqRaw.vcf.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(RAW,stri_locate_last_fixed(RAW, "/")[,1]+1, nchar(RAW))) %>%
  mutate(SampleID_Flowcell = sub(".lofreqRaw.vcf.gz", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)


SEXs <- data.frame(SEX = list.files(paste0(nxf_work_dir, "/../QC_",version,"/"), "*.predictedSex.txt$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(SEX,stri_locate_last_fixed(SEX, "/")[,1]+1, nchar(SEX))) %>%
  mutate(SampleID_Flowcell = sub(".predictedSex.txt", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)

mtDNAs <- data.frame(mtDNA = list.files(paste0(nxf_work_dir, "/../QC_",version,"/"), "*.mtDNAhg_classified.txt$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(mtDNA,stri_locate_last_fixed(mtDNA, "/")[,1]+1, nchar(mtDNA))) %>%
  mutate(SampleID_Flowcell = sub(".mtDNAhg_classified.txt", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)

ANCESTRYs <- data.frame(ANCESTRY = list.files(paste0(nxf_work_dir, "/../QC_",version,"/"), "*.ancestryPrediction.txt$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(ANCESTRY,stri_locate_last_fixed(ANCESTRY, "/")[,1]+1, nchar(ANCESTRY))) %>%
  mutate(SampleID_Flowcell = sub(".ancestryPrediction.txt", "", SampleID_Flowcell)) %>% filter(SampleID_Flowcell %in% NBSNGS_sampleSheet$SampleID_Flowcell)


#Creating connection: ----
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/variants.sqlite"))
#KILL
    #dbExecute(con, "DROP TABLE IF EXISTS LOFREQ;")
    #dbExecute(con, "DROP TABLE IF EXISTS OPL;")
    #dbExecute(con, "DROP TABLE IF EXISTS COV;")
    #dbExecute(con, "DROP TABLE IF EXISTS RAW;")
    #dbExecute(con, "DROP TABLE IF EXISTS SEX;")
    #dbExecute(con, "DROP TABLE IF EXISTS mtDNA;")
    #dbExecute(con, "DROP TABLE IF EXISTS ANCESTRY;")


#BUILD database: ----
#Check if databases exists:
LOFREQ_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='LOFREQ';")
OPL_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='OPL';")
COV_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='COV';")
RAW_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='RAW';")
SEX_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='SEX';")
mtDNA_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='mtDNA';")
ANCESTRY_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='ANCESTRY';")

#OPL SQL ----
if (nrow(OPL_exists) > 0) {
existing_opl_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM OPL")
add_OPL <- OPLs %>% filter(!SampleID_Flowcell %in% existing_opl_SampleID_Flowcells$SampleID_Flowcell)
if(nrow(add_OPL) > 0){
  print(paste0("OPL: adding content from ", nrow(add_OPL), " samples!"))
}else{
  print("OPL is already up to date!")
}  
}else{
  add_OPL <- OPLs
  print(paste0("OPL did NOT exists - bulding database from ", nrow(add_OPL), " samples!"))
  dbExecute(con, "
  CREATE TABLE OPL (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT ,
    CHROM TEXT NOT NULL,
    POS INTEGER ,
    REF TEXT NOT NULL ,
    ALT TEXT NOT NULL ,
    RSID TEXT,
    AF INTEGER,
    FS INTEGER, 
    AF_nfe INTEGER, 
    CLNSIG TEXT,
    RS TEXT, 
    CLNVC TEXT, 
    CLNDN TEXT, 
    CLNREVSTAT TEXT,
    GENEINFO TEXT,
    GENEID TEXT,
    HGVS_P TEXT, 
    FEATUREID TEXT, 
    EFFECT TEXT, 
    GENE TEXT,
    HGVS_C TEXT,
    GT TEXT, 
    FWD TEXT, 
    REV TEXT,
    DP INTEGER,
    GQ INTEGER, 
    PL TEXT
     )
") } 
if(nrow(add_OPL) > 0){
for(i in 1:nrow(add_OPL)){
  OPL <- fread(cmd = paste0("zcat ", add_OPL$OPL[i]), na.strings = "") %>% 
    mutate(SampleID_Flowcell = add_OPL$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything()) %>% rename(RSID = ID) %>% 
    separate(`GEN[*]`, c("GT", "AD", "DP", "GQ", "PL"), sep = ":") %>% separate(AD, c("FWD","REV"), sep = ",") %>% 
    mutate(AF = as.numeric(REV)/(as.numeric(FWD)+as.numeric(REV))) %>% mutate(across(c(GQ, DP), as.numeric))  %>% mutate(AF_nfe = as.numeric(AF_nfe)) %>% replace_na(list(AF_nfe = 0))
  colnames(OPL) <- gsub("ANN\\[\\*\\]\\.", "", colnames(OPL))
  OPL_wrangled <- OPL %>% group_by(SampleID_Flowcell, CHROM, POS, REF, ALT) %>% summarise(
    across(.cols  = everything(),.fns   = ~ paste(unique(na.omit(.)), collapse = ":"),.names = "{.col}" ),.groups = "drop")
  if(nrow(OPL_wrangled) == 0){
    OPL_wrangled <- tibble(!!!setNames(rep(list(NA), length(OPL)), names(OPL))) %>% mutate(SampleID_Flowcell = add_OPL$SampleID_Flowcell[i], CHROM = "NA", REF = "NA", ALT = "NA")
  }
  dbAppendTable(con, "OPL", OPL_wrangled, conflict = "replace")
  }
}  

#LOFREQ SQL ----
if (nrow(LOFREQ_exists) > 0) {
  existing_lofreq_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM LOFREQ")
  add_LOFREQ <- LOFREQs %>% filter(!SampleID_Flowcell %in% existing_lofreq_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_LOFREQ) > 0){
  print(paste0("LOFREQ: adding content from ", nrow(add_LOFREQ), " samples!"))
  }else{
    print("LOFREQ is already up to date!")
  }  
}else{
  add_LOFREQ <- LOFREQs
  print(paste0("LOFREQ did NOT exists - bulding database from ", nrow(add_LOFREQ), " samples!"))
dbExecute(con, "
  CREATE TABLE LOFREQ (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    VARID TEXT NOT NULL,
    POS INTEGER,
    AF INTEGER,
    DP INTEGER, 
    SB INTEGER,
    Fref TEXT, 
    Rref TEXT,
    Falt TEXT,
    Ralt TEXT
     )
")}
if(nrow(add_LOFREQ) > 0){
for(i in 1:nrow(add_LOFREQ)){
  LOFREQ  <- fread(cmd=paste0("zcat ",add_LOFREQ$LOFREQ[i],"|  grep -v '^##'" ))  %>% rename(CHROM = 1) %>% 
    unite(VARID, c(CHROM,POS, REF, ALT), sep = "-", remove = FALSE) %>%
    separate(., INFO, c("DP", "AF", "SB", "DP4"), sep = ";")  %>% 
    mutate(DP = gsub("DP=", "", DP), AF = gsub("AF=", "" , AF), SB = gsub("SB=", "" , SB), DP4 = gsub("DP4=", "" , DP4) ) %>% 
    separate(DP4, c("Fref","Rref", "Falt", "Ralt"), sep = ",") %>% mutate(across(c(AF, SB, DP), as.numeric)) %>% 
    select(VARID, POS, AF, DP, SB, Fref, Rref, Falt, Ralt)  %>%  
    mutate(SampleID_Flowcell = add_LOFREQ$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything())
  if(nrow(LOFREQ) == 0){
    LOFREQ <- tibble(!!!setNames(rep(list(NA), length(LOFREQ)), names(LOFREQ))) %>% mutate(SampleID_Flowcell = add_LOFREQ$SampleID_Flowcell[i], VARID = "NA")
  }
  dbAppendTable(con, "LOFREQ", LOFREQ, conflict = "replace")
  }
}
#COV SQL ----
if (nrow(COV_exists) > 0) {
  existing_cov_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM COV")
  add_COV <- COVs %>% filter(!SampleID_Flowcell %in% existing_cov_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_COV) > 0){
    print(paste0("COV: adding content from ", nrow(add_COV), " samples!"))
  }else{
    print("COV is already up to date!")
  }  
}else{
  add_COV <- COVs
  print(paste0("COV did NOT exists - bulding database from ", nrow(add_COV), " samples!"))
  dbExecute(con, "
  CREATE TABLE COV (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    CHR_POS TEXT,
    DP INTEGER 
     )
")}
if(nrow(add_COV) > 0){
  for(i in 1:nrow(add_COV)){
    COVtest <- fread(cmd = paste0("zcat ", add_COV$COV[i], "| wc -l")) %>% pull(V1)
    if(COVtest > 1){
    COV <- fread(cmd = paste0("zcat ", add_COV$COV[i])) %>% unite(CHR_POS, 1:2, sep = "-") %>% rename(DP = 2) %>% 
      filter(CHR_POS %in% bedfile_coverage_expanded_targets$chr.pos) %>% 
      mutate(SampleID_Flowcell = add_COV$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything())
    }else{
     COV <- bedfile_coverage_expanded_targets  %>% select(chr.pos) %>% mutate(DP = 0) %>% rename(CHR_POS = chr.pos) %>% mutate(SampleID_Flowcell = add_COV$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything())
    }
    dbAppendTable(con, "COV", COV, conflict = "replace")
    }
}  

#RAW SQL ----
if (nrow(RAW_exists) > 0) {
  existing_raw_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM RAW")
  add_RAW <- RAWs %>% filter(!SampleID_Flowcell %in% existing_raw_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_RAW) > 0){
    print(paste0("RAW: adding content from ", nrow(add_RAW), " samples!"))
  }else{
    print("RAW is already up to date!")
  }  
}else{
  add_RAW <- RAWs
  print(paste0("RAW did NOT exists - bulding database from ", nrow(add_RAW), " samples!"))
  dbExecute(con, "
  CREATE TABLE RAW (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    VARID TEXT,
    SB TEXT,
    DP4 TEXT
     )
")}

if(nrow(add_RAW) > 0){
for(i in 1:nrow(add_RAW)){
RAW <- fread(cmd = paste0(
  "zcat ", add_RAW$RAW[i], 
  " | grep -v '^##' | awk 'NR==1 {print $0} NR>1 && $8 ~ /DP4=/ && $8 ~ /AF=/ {match($8, /DP4=([0-9]+),([0-9]+),([0-9]+),([0-9]+)/, dp4); z=dp4[3]; a=dp4[4]; match($8, /AF=([0-9\\.]+)/, af); if (z+a >= 3 && af[1] >= 0.025) print $0}'")) %>% 
  rename(CHROM = 1) %>% unite(VARID, c(CHROM,POS, REF, ALT), sep = "-") %>%
  separate(., INFO, c("DP", "AF", "SB", "DP4"), sep = ";")  %>% 
  mutate(DP = gsub("DP=", "", DP), AF = gsub("AF=", "" , AF), SB = gsub("SB=", "" , SB), DP4 = gsub("DP4=", "" , DP4) ) %>% 
  select(VARID, SB, DP4) %>% mutate(SampleID_Flowcell = add_RAW$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything())
if(nrow(RAW) == 0){
  RAW <- tibble(!!!setNames(rep(list(NA), length(RAW)), names(RAW))) %>% mutate(SampleID_Flowcell = add_RAW$SampleID_Flowcell[i])
}
dbAppendTable(con, "RAW", RAW, conflict = "replace")

  }
}

#Sex table
if (nrow(SEX_exists) > 0) {
  existing_sex_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM SEX")
  add_SEX <- SEXs %>% filter(!SampleID_Flowcell %in% existing_sex_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_SEX) > 0){
    print(paste0("SEX: adding content from ", nrow(add_SEX), " samples!"))
  }else{
    print("SEX is already up to date!")
  }  
}else{
  add_SEX <- SEXs
  print(paste0("SEX did NOT exists - bulding database from ", nrow(add_SEX), " samples!"))
  dbExecute(con, "
  CREATE TABLE SEX (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    cov_chrX TEXT,
    cov_chrY TEXT,
    XY_ratio TEXT,
    Predicted_Sex TEXT
     )
")}
if(nrow(add_SEX) > 0){
  for(i in 1:nrow(add_SEX)){
sex <- fread(add_SEX$SEX[i]) %>%rename(cov_chrX = chrX, cov_chrY = chrY, XY_ratio = XY.ratio, Predicted_Sex = Sex.prediction) %>% 
  select(SampleID_Flowcell,cov_chrX , cov_chrY,XY_ratio,Predicted_Sex)
dbAppendTable(con, "SEX", sex, conflict = "replace")
  }
}



#mtDNA table
if (nrow(mtDNA_exists) > 0) {
  existing_mtDNA_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM mtDNA")
  add_mtDNA <- mtDNAs %>% filter(!SampleID_Flowcell %in% existing_mtDNA_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_mtDNA) > 0){
    print(paste0("mtDNA: adding content from ", nrow(add_mtDNA), " samples!"))
  }else{
    print("mtDNA is already up to date!")
  }  
}else{
  add_mtDNA <- mtDNAs
  print(paste0("mtDNA did NOT exists - bulding database from ", nrow(add_mtDNA), " samples!"))
  dbExecute(con, "
  CREATE TABLE mtDNA (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    Haplogroup TEXT,
    Quality TEXT,
    Variants TEXT
     )
")}
if(nrow(add_mtDNA) > 0){
  for(i in 1:nrow(add_mtDNA)){
    MTDNA <- fread(add_mtDNA$mtDNA[i]) %>% select(SampleID, Haplogroup, Quality, length(.)) %>% rename(SampleID_Flowcell = SampleID, Variants = 4)
    dbAppendTable(con, "mtDNA", MTDNA, conflict = "replace")
  }
}



#ANCESTRY table
if (nrow(ANCESTRY_exists) > 0) {
  existing_ANCESTRY_SampleID_Flowcells <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM ANCESTRY")
  add_ANCESTRY <- ANCESTRYs %>% filter(!SampleID_Flowcell %in% existing_ANCESTRY_SampleID_Flowcells$SampleID_Flowcell)
  if(nrow(add_ANCESTRY) > 0){
    print(paste0("ANCESTRY: adding content from ", nrow(add_ANCESTRY), " samples!"))
  }else{
    print("ANCESTRY is already up to date!")
  }  
}else{
  add_ANCESTRY <- ANCESTRYs
  print(paste0("ANCESTRY did NOT exists - bulding database from ", nrow(add_ANCESTRY), " samples!"))
  dbExecute(con, "
  CREATE TABLE ANCESTRY (
    id INTEGER PRIMARY KEY,
    SampleID_Flowcell TEXT,
    INFO TEXT
     )
")}

if(nrow(add_ANCESTRY) > 0){
  for(i in 1:nrow(add_ANCESTRY)){
    ANCESTRY_import <- fread(add_ANCESTRY$ANCESTRY[i]) %>% mutate(Accuracy = round(as.numeric(Accuracy),digits = 2)) %>% 
      mutate(INFO = paste0(Ancestry," (Probability=", Probability,";Info=", Info, ";Accuracy=", Accuracy,")")) %>% select(-Probability, -Info, -Accuracy, -Ancestry)
    dbAppendTable(con, "ANCESTRY", ANCESTRY_import, conflict = "replace")
  }
}



dbDisconnect(con)


OPLlist <- OPLs$SampleID_Flowcell %>% paste(., collapse = "|")
done.REPORTS <-  data.frame(REPORTS =  list.files(paste0(nxf_work_dir, "/../REPORTS_",version,"/"), "*.report.html$", recursive = TRUE)) %>% 
  mutate(SampleId_Flowcell =  str_extract(REPORTS, OPLlist  ))

samplesForProcessing <- NBSNGS_sampleSheet %>% filter(!SampleID_Flowcell %in% done.REPORTS$SampleId_Flowcell ) %>% 
  select(SampleID_Flowcell, Flowcell, Description) 


fwrite(samplesForProcessing, "SAMPLESforREPORTING.txt", sep = "\t")



#