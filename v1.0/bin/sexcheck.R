#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
version <- args[2]

#nxf_work_dir <-  "~/TEST_TMP/NBSNGS_test/NBSNGS_WF/v1.0"
#version <- "v.1.0"

source(paste0(nxf_work_dir,"/bin/dependencies.R"))


bed_expanded_coverage <- fread( paste0(nxf_work_dir,"/assets/temp/bedfile_coverage_expanded.bed")) %>% filter(grepl("chrX|chrY",chr.pos)) %>% separate(chr.pos, c("Chromosome", "pos"), sep ="-") %>% 
  mutate(pos = as.integer(pos))

COVERAGE <- data.frame(COVERAGE = list.files(paste0(nxf_work_dir,"/../QC_",version,"/"), "*.cov.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(COVERAGE,stri_locate_last_fixed(COVERAGE, "/")[,1]+1, nchar(COVERAGE))) %>%
  mutate(SampleID_Flowcell  = sub(".cov.gz", "", SampleID_Flowcell ))  %>% separate(SampleID_Flowcell, c("SampleID", "Flowcell"), sep = "-", remove = FALSE)
existing_sexcheck <- data.frame(Sex = list.files(paste0(nxf_work_dir,"/../QC_",version,"/"),"*.predictedSex.txt", recursive = TRUE, full.names = TRUE)) %>% 
  mutate(SampleID_Flowcell = str_sub(Sex,stri_locate_last_fixed(Sex, "/")[,1]+1, nchar(Sex))) %>%
  mutate(SampleID_Flowcell  = sub(".predictedSex.txt", "", SampleID_Flowcell ))  

DB <-  COVERAGE %>% filter(!SampleID_Flowcell %in% existing_sexcheck$SampleID_Flowcell) 
#i <- 2
#collect_results <- list()
if(nrow(DB) >= 1){
for(i in 1:nrow(DB)){
  print(i)
  importCov <- fread(cmd = paste0("zcat ", DB$COVERAGE[i], " | grep chrX | wc -l"))  %>% pull(V1)
    if(importCov > 1){
    chrXY <- fread(DB$COVERAGE[i]) %>% rename(Chromosome = 1, pos = 2)%>% inner_join(., bed_expanded_coverage, by = c("Chromosome", "pos")) %>% 
    filter(Chromosome %in% c("chrX", "chrY")) %>% group_by(Chromosome) %>%  summarise(Mean = mean(V3)) %>% spread(Chromosome, Mean)%>% 
    mutate(chrX= ifelse(!("chrX" %in% names(.)), 0, ifelse(is.na(chrX), 0, chrX)),
    chrY = ifelse(!("chrY" %in% names(.)), 0, ifelse(is.na(chrY), 0, chrY))) %>% mutate( SampleID_Flowcell = DB$SampleID_Flowcell[i]) %>% 
    mutate(XY.ratio = chrX/chrY) %>% mutate(Sex.prediction = ifelse(XY.ratio < 3, "Male", "Female")) %>% 
    mutate(Sex.prediction = ifelse(is.na(Sex.prediction), "Insufficent data", Sex.prediction))
    }else{
    chrXY <- data.frame(chrX = 0, chrY = 0)%>% mutate( SampleID_Flowcell = DB$SampleID_Flowcell[i]) %>% 
    mutate(XY.ratio = chrX/chrY) %>% mutate(Sex.prediction = ifelse(XY.ratio < 3, "Male", "Female")) %>% 
    mutate(Sex.prediction = ifelse(is.na(Sex.prediction), "Insufficent data", Sex.prediction))
    }
  fwrite( chrXY, paste0(nxf_work_dir,"/../QC_",version,"/",DB$Flowcell[i],"/",DB$SampleID_Flowcell[i],".predictedSex.txt"), sep = "\t")
  }
}

if(nrow(DB) == 0){
DB <- data.frame(VCF = "-", SampleID_Flowcell = "-", SampleID = "-", Flowcell = "-", COVERAGE = "-")
}
fwrite(DB, "SexProcessedIDs.txt", sep = "\t")

