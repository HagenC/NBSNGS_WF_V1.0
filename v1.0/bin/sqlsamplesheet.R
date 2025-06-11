#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
bedfile_coverage <- args[2]
bedfile_call <- args[3]
bedfile_target <- args[4]
snpeff_gtf <- args[5]
version <- args[6]

#nxf_work_dir <-  "~/TEST_TMP/NBSNGS_test/NBSNGS_WF/v1.0/"
#version <- "v.1.0"

source(paste0(nxf_work_dir,"/bin/dependencies.R"))

#List files  ----
sample.sheet.list <- list.files(paste0(nxf_work_dir), "^SampleSheet.csv$", recursive = TRUE, full.names = TRUE)
copy.complete.list <- list.files(paste0(nxf_work_dir ), "CopyComplete.txt$", recursive = TRUE, full.names = T) %>% gsub("CopyComplete.txt", "SampleSheet.csv", .) 
copy.complete.list <- grep("NBSNGS.PIPELINE|arkiv", copy.complete.list, invert = TRUE, value = TRUE)
sample.sheet.processing.list <- sample.sheet.list[sample.sheet.list %in% copy.complete.list]
#sample.sheet.list -> sample.sheet.processing.list

#Importing SampleSheets ----
position <- 1
sample.sheet.collect <- list()
for(i in 1:length(sample.sheet.processing.list)){
  names_check <- names(fread(cmd=paste0("sed -n '/Sample_ID/,$p' ", sample.sheet.processing.list[i]), nrows = 1))
  required_names <- c("Sample_ID", "Sample_Name", "Sample_Project", "Description")
  if (all(required_names %in% names_check)) {
  sample.sheet.collect[[position]] <- fread(cmd=paste0("sed -n '/Sample_ID/,$p' ",  sample.sheet.processing.list[i]), na.strings = "") %>% 
    select(Sample_ID,Sample_Name, Sample_Project, Description) %>% mutate(SampleSheet.path = sample.sheet.processing.list[i]) %>% 
    mutate(Flowcell = str_extract(SampleSheet.path, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+")) %>% filter(!is.na(Flowcell)) %>% 
    unite(Sample_ID.Flowcell, c(Sample_ID, Flowcell), sep = "-", remove= FALSE)
  position <- position +1
  }
}
sampleSheets <- do.call("rbind", sample.sheet.collect) %>% filter(!is.na(Sample_ID)) %>% distinct() %>% filter(!is.na(Sample_Project))

#Searching for fastq.gz files ----
fastQ.list <- sample.sheet.processing.list %>% gsub("SampleSheet.csv", "",.) 
sample.fastq.collect <- list()
for(i in 1:length(fastQ.list)){
  sample.fastq.collect[[i]] <-  data.frame(FQs = list.files( fastQ.list[i], "*.fastq.gz$", recursive = TRUE, full.names = T))
}


fastqs <- do.call("rbind", sample.fastq.collect) %>%
  mutate(File = str_sub(FQs,stri_locate_last_fixed(FQs, "/")[,1]+1, nchar(FQs))) %>% 
  mutate(Sample_ID = gsub("_[[:alpha:]][[:digit:]]+_[[:alpha:]][[:digit:]]_001\\.fastq\\.gz","",File)) %>%  rowwise() %>%  
  mutate(flowcell = str_extract(FQs, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+")) %>%  mutate(Sample_ID.Flowcell = paste0(Sample_ID, "-",flowcell)) %>% 
  filter(Sample_ID.Flowcell %in% sampleSheets$Sample_ID.Flowcell) %>% select(-File)

sample.list.sizeInfo <- data.frame(file.info(fastqs$FQs)) 

#Fetching size and mtime infor for fastq.gz files ----
sample.size <- sample.list.sizeInfo %>% mutate(fastq.size.gb = round((size/(1024^2)/1000), digits = 3)) %>% select(fastq.size.gb) %>% 
  mutate(sample_path = row.names(.)) %>% mutate(File = str_sub(sample_path,stri_locate_last_fixed(sample_path, "/")[,1]+1, nchar(sample_path))) %>% 
  mutate(Sample_ID = gsub("_[[:alpha:]][[:digit:]]+_[[:alpha:]][[:digit:]]_001\\.fastq\\.gz","",File)) %>% 
  mutate(flowcell = str_extract(sample_path, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+"))  %>%
  mutate(Sample_ID.Flowcell = paste0(Sample_ID, "-",flowcell))%>% select(Sample_ID.Flowcell, fastq.size.gb)%>% filter(Sample_ID.Flowcell %in% sampleSheets$Sample_ID.Flowcell) %>%  
  group_by(Sample_ID.Flowcell) %>%  mutate(order = seq_along(Sample_ID.Flowcell)) %>% spread(order, fastq.size.gb) %>% 
  rename(fastq.1.size.GB = 2, fastq.2.size.GB = 3)
sample.time <-sample.list.sizeInfo %>%  select(mtime) %>% mutate(sample_path = row.names(.)) %>% 
   mutate(File = str_sub(sample_path,stri_locate_last_fixed(sample_path, "/")[,1]+1, nchar(sample_path))) %>% 
  mutate(Sample_ID = gsub("_[[:alpha:]][[:digit:]]+_[[:alpha:]][[:digit:]]_001\\.fastq\\.gz","",File)) %>% 
  mutate(flowcell = str_extract(sample_path, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+"))  %>%
  mutate(Sample_ID.Flowcell = paste0(Sample_ID, "-",flowcell)) %>%  select(Sample_ID.Flowcell, mtime)%>% filter(Sample_ID.Flowcell %in% sampleSheets$Sample_ID.Flowcell) %>%  
  group_by(Sample_ID.Flowcell) %>%  mutate(order = seq_along(Sample_ID.Flowcell)) %>% spread(order, mtime) %>% 
  rename(mtime = 2, mtime.2 = 3) %>% select(-mtime.2) 

sample.list.sizeInfo.time <- left_join(sample.size, sample.time, by = "Sample_ID.Flowcell")

# Combining SampleSheet, Fastq files and info ----
#Filtering on non missing fastq.gz. files. 
SampleSheet_fastq_DB <- fastqs %>%  rowwise() %>% 
  ungroup() %>% group_by(Sample_ID.Flowcell) %>%  mutate(order = seq_along(Sample_ID.Flowcell)) %>% spread(order, FQs) %>% select(-Sample_ID) %>% 
  right_join(., sampleSheets, by = c("Sample_ID.Flowcell")) %>% 
  left_join(.,sample.list.sizeInfo.time, by = "Sample_ID.Flowcell" ) %>% 
  arrange(mtime) %>% 
  rename(fastq_1 = `1`, fastq_2 = `2`,samplesheet_path = SampleSheet.path, fastq_1_sizeGB = fastq.1.size.GB,  fastq_2_sizeGB = fastq.2.size.GB, SampleID_Flowcell = Sample_ID.Flowcell) %>% 
  select(Sample_ID, Flowcell, SampleID_Flowcell, Description, fastq_1, fastq_2, everything(), - flowcell) %>% filter(!is.na(fastq_2) & !is.na(fastq_1)) 



#Building SQL database: ----
#Creating connection: 
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir, "/assets/sql/SampleSheets.sqlite"))
#KILL
  #dbExecute(con, "DROP TABLE IF EXISTS SAMPLESHEETS;")

#Check if databases exists:
SAMPLESHEETS_exists <- dbGetQuery(con,"SELECT name FROM sqlite_master WHERE type='table' AND name='SAMPLESHEETS';")

if (nrow(SAMPLESHEETS_exists) > 0) {
  existing_samplesheets_SampleID_Flowcell <- dbGetQuery(con, "SELECT DISTINCT SampleID_Flowcell FROM SAMPLESHEETS")
  add_SAMPLESHEETS <- SampleSheet_fastq_DB  %>% filter(!SampleID_Flowcell %in% existing_samplesheets_SampleID_Flowcell$SampleID_Flowcell)
  print(paste0("SAMPLESHEETS exists - adding content from ", nrow(add_SAMPLESHEETS), " samples!"))
}else{
  add_SAMPLESHEETS  <- SampleSheet_fastq_DB
  print(paste0("SAMPLESHEETS did NOT exists - bulding database from ", nrow(add_SAMPLESHEETS), " samples!"))
  dbExecute(con, "
            CREATE TABLE SAMPLESHEETS (
            id INTEGER PRIMARY KEY,
            Sample_ID TEXT, 
            Flowcell TEXT ,
            SampleID_Flowcell TEXT NOT NULL UNIQUE,
            Description TEXT  ,
            fastq_1 TEXT NOT NULL ,
            fastq_2 TEXT NOT NULL ,
            Sample_Name TEXT,
            Sample_Project TEXT,
            samplesheet_path TEXT, 
            fastq_1_sizeGB REAL, 
            fastq_2_sizeGB REAL, 
            mtime TIMESTAMP
                    )
            ") 
  } 
if (nrow(add_SAMPLESHEETS) > 0) {
dbAppendTable(con, "SAMPLESHEETS", add_SAMPLESHEETS, conflict = "replace")
}
dbDisconnect(con)  


con <- dbConnect(RSQLite::SQLite(),  dbname = paste0(nxf_work_dir, "/assets/sql/SampleSheets.sqlite"))
sampleProject <- c("NBS-NGS", "NBS_NGS", "nbs-ngs") %>% paste(collapse = "', '") %>%   paste0("'", ., "'") 
 query_samplesheets <- paste0("
   SELECT *
   FROM SAMPLESHEETS
   WHERE Sample_Project IN (",sampleProject,")
 ")
 
NBSNGS_samples <- dbGetQuery(con, query_samplesheets)
dbDisconnect(con) 

#Done:
DONE_gvcfs <- data.frame(GVCF =  list.files(paste0(nxf_work_dir, "/../variants_",version,"/"), "*.GATK.g.vcf.gz$", recursive = TRUE)) %>% 
  mutate(Sample_ID.gvcf = str_sub(GVCF,stri_locate_last_fixed(GVCF, "/")[,1]+1, nchar(GVCF))) %>% mutate(SampleID_Flowcell = gsub("*.GATK.g.vcf.gz", "", Sample_ID.gvcf))


#Alignemt tuple 
N = 11
max <- nrow(NBSNGS_samples)
samplesTorun <- ifelse(max >= N, N, max)
alignment_list <- NBSNGS_samples %>% arrange(desc(mtime)) %>%  filter(!grepl("JOGR|HCP", Sample_ID)) %>% select(SampleID_Flowcell,  fastq_1, fastq_2, Flowcell) %>% filter(!SampleID_Flowcell %in% DONE_gvcfs$SampleID_Flowcell)
if(nrow(alignment_list) > N){
  alignment_list <- alignment_list %>% .[1:samplesTorun,]
}

#Mail attachment
mailAttach <- NBSNGS_samples %>% filter(SampleID_Flowcell %in% alignment_list$SampleID_Flowcell) %>% select(Sample_ID, Sample_Name, Description, SampleID_Flowcell, 
                                                                                                                    fastq_1_sizeGB, fastq_2_sizeGB, mtime)

#Export
fwrite(NBSNGS_samples, "Sample.collection.info", sep = "\t")
fwrite(mailAttach, "ProcessedSamples.txt", sep = "\t")
fwrite(alignment_list, "ProccessingIDs.txt", sep = "\t")

#asseets prep ----

#SnpEff gene filter
phenotype_list <- fread(paste0(nxf_work_dir,"/assets/phenotype/phenotype.gene.list")) %>% .[,2]
fwrite(phenotype_list, paste0(nxf_work_dir,"/assets/phenotype/geneFilter.list"), col.names = FALSE)

foo <- fread("/srv/data/bedFiles/NBSNGS_v2/NBSNGS_v2_CALL.bed")
bedfile_target <- "/srv/data/bedFiles/NBSNGS_v2/NBSNGS_v2_CALL.bed"

#Expanded beds:
if(file.exists(paste0(nxf_work_dir,"/assets/temp/bedfile_call_expanded.bed"))== "FALSE"){
  bedfile_coverage_i <- fread(bedfile_coverage) %>% rename(Chromosome = 1, Start = 2, End = 3,GeneID = 5, Exon = 4) %>% select(1:5) 
  bedfile_call_i <-  fread(bedfile_call) %>% rename(Chromosome = 1, Start = 2, End = 3,GeneID = 5, Exon = 4) %>% select(1:5) %>% mutate(Size = End - Start)
  bedfile_targets_i <-  fread(bedfile_target) %>% rename(Chromosome = 1, Start = 2, End = 3,GeneID = 5, Exon = 4) %>% select(1:5)  %>% 
    filter(GeneID %in% phenotype_list$Gene.ID)
  expandFunction <- function(N) {
    expandedRegion <- N$Start:N$End
    return(data.frame(chr = N$Chromosome,  pos = expandedRegion, Start = N$Start, End = N$End, GeneID = N$GeneID, Exon = N$Exon)) #, Size = N$End - N$Start
  }
  
  bedfile_coverage_expanded <- do.call(rbind, lapply(1:nrow(bedfile_coverage_i), function(i) expandFunction(bedfile_coverage_i[i, ]))) %>% unite(chr.pos, chr:pos, sep = "-") 
  bedfile_call_expanded <- do.call(rbind, lapply(1:nrow(bedfile_call_i), function(i) expandFunction(bedfile_call_i[i, ]))) %>% unite(chr.pos, chr:pos, sep = "-") %>% filter(!duplicated(chr.pos))
  bedfile_targets <- do.call(rbind, lapply(1:nrow(bedfile_targets_i), function(i) expandFunction(bedfile_targets_i[i, ]))) %>% unite(chr.pos, chr:pos, sep = "-") %>% filter(!duplicated(chr.pos))
  fwrite(bedfile_targets , paste0(nxf_work_dir,"/assets/temp/bedfile_target.bed"), sep = "\t")
  fwrite(bedfile_coverage_expanded , paste0(nxf_work_dir,"/assets/temp/bedfile_coverage_expanded.bed"), sep = "\t")
  fwrite(bedfile_call_expanded, paste0(nxf_work_dir,"/assets/temp/bedfile_call_expanded.bed"), sep = "\t")
}





#Sanity check of SnpEff.gtf database i.e /assets/phenotype/geneFilter.list is used to filter taget genes and needs to match SnpEff annotations. 
if(file.exists(paste0(nxf_work_dir,"/assets/sanity/SnpEff_geneFilter.list_observation_count.tsv"))== "FALSE"){
collect_SnpEff_database_gtf_counts  <- list()
for (i in 1:nrow(phenotype_list)) {
  gene <- phenotype_list$Gene.ID[i]
  pattern <- sprintf('gene_name "%s"', gene)
  collect_SnpEff_database_gtf_counts[[i]] <- data.frame(Gene.ID = phenotype_list$Gene.ID[i], SnpEff.observations.gtf = system(paste0("grep -F -c ",shQuote(pattern)," ",snpeff_gtf), intern = TRUE))
}

SnpEff_observations_counts <- do.call("rbind",collect_SnpEff_database_gtf_counts)
fwrite(SnpEff_observations_counts, paste0(nxf_work_dir,"/assets/sanity/SnpEff_geneFilter.list_observation_count.tsv"), sep = "\t")
}
