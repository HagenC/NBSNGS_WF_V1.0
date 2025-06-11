#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]


#nxf_work_dir <-  "/srv/data/ILMN/PIPELINES/NbsNgsWF"
source(paste0(nxf_work_dir,"/bin/dependencies.R"))

LA <- data.frame(LA = list.files("/srv/data/ILMN/PIPELINES/NbsNgsWF/variants/", "*.GATK.LA.vcf.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(LA,stri_locate_last_fixed(LA, "/")[,1]+1, nchar(LA))) %>%
  mutate(SampleID_Flowcell  = sub(".GATK.LA.vcf.gz", "", SampleID_Flowcell )) %>%
  mutate(Flowcell = str_extract(SampleID_Flowcell, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+")) %>% select(SampleID_Flowcell,Flowcell,LA) # %>% .[1:10,]

fwrite(LA, "LA.tsv", sep =  "\t")






