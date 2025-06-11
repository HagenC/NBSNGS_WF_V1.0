#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
ancestry_reference_data <- args[2]
ancestry_reference_RS <- args[3]
version <- args[4]

#nxf_work_dir <-  "/srv/data/ILMN/PIPELINES/NBSNGS_WF_v1.0/v1.0/"
#version <- "v.1.0"
#setwd(nxf_work_dir)
#ancestry_reference_RS <- "/srv/data/VCF_annotations/Ancestry_prediction/reference_RS.tsv"
#ancestry_reference_data <- "/srv/data/VCF_annotations/Ancestry_prediction/Ancestry_reference_data.tsv"
source(paste0(nxf_work_dir,"/bin/dependencies.R"))


#Importing sampleheets: ----
con <- dbConnect(RSQLite::SQLite(), dbname = paste0(nxf_work_dir , "/assets/sql/SampleSheets.sqlite"))
sampleProject <- c("NBS-NGS", "nbs-ngs", "NBS_NGS") %>% paste(collapse = "', '") %>%   paste0("'", ., "'") 
query_opl <- paste0(" SELECT * FROM SAMPLESHEETS WHERE Sample_Project IN (",sampleProject,") ")
NBSNGS_sampleSheet <- dbGetQuery(con, query_opl) %>% select(SampleID_Flowcell, Description)
dbDisconnect(con) 

DONE <- data.frame(predicted = list.files(paste0(nxf_work_dir,"/../QC_",version,"/"), "*.ancestryPrediction.txt", recursive = TRUE)) %>%
  mutate(predicted = gsub(".ancestryPrediction.txt", "", predicted)) %>% 
  mutate(SampleID_Flowcell = str_sub(predicted,stri_locate_last_fixed(predicted, "/")[,1]+1, nchar(predicted)))
VCF_list <- data.frame(VCF = list.files(paste0(nxf_work_dir,"/../variants_",version,"/"), "*GATK.g.vcf.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(VCF,stri_locate_last_fixed(VCF, "/")[,1]+1, nchar(VCF))) %>%
  mutate(SampleID_Flowcell  = sub(".GATK.g.vcf.gz", "", SampleID_Flowcell )) %>% 
  separate(SampleID_Flowcell, c("SampleID", "Flowcell"), sep = "-", remove = FALSE) %>% 
  filter(!SampleID_Flowcell %in% DONE$SampleID_Flowcell) 

#i <- 1
set.seed(42)
if(nrow(VCF_list) > 0){
RS_selection <- fread(ancestry_reference_RS)
referenceData <- fread(ancestry_reference_data) %>% mutate(Ancestry = as.factor(Ancestry))
for(i in 1:nrow(VCF_list)){
  sampleVCF <- fread(cmd = paste0("zcat ",VCF_list$VCF[i],"  | grep  -v \"^##\" | awk '$8 !~ /END=[0-9]+/'")) %>% mutate(`#CHROM` = gsub("chr", "", `#CHROM`), ALT = gsub(",<NON_REF>","", ALT)) %>%  
    rename(Sample = 10) %>% 
    mutate(`#CHROM` = as.integer(`#CHROM`)) %>% filter(!is.na(`#CHROM`)) %>% mutate(VARID = paste0(`#CHROM`,"-", POS,"-", REF,"-", ALT )) 
  INFO <- sampleVCF %>% filter(VARID %in% RS_selection$VARID)
  if(nrow(INFO) > 100){
  conversionTable <- data.frame(GT = c("0/0", "0/1","1/1", "1|1","0|1", "0|0", "1|0"), value = c(0,1,2,2,1,0,1))
  samplePredVCF  <- inner_join(sampleVCF, RS_selection, by = "VARID") %>%  
    mutate(ID = RS) %>% 
    select(ID, Sample)  %>% mutate(GT = str_extract(Sample, "0/0|0/1|1/1|1\\|1|0\\|1")) %>% filter(!is.na(GT)) %>% select(ID, GT) %>% left_join(., conversionTable, by = "GT") %>%
    select(-GT) %>% mutate(Ancestry = VCF_list$SampleID_Flowcell[i]) %>% mutate(ID = gsub(";", ".", ID)) %>%  pivot_wider(
      names_from = ID,  values_from = value,values_fill = NA) 
  VARIANTS <- length(samplePredVCF) - 1
  modelData <- referenceData %>% select(colnames(samplePredVCF))
  #Split for Accuracy measure
  sample_indices <- sample(1:nrow(modelData), 0.8 * nrow(modelData))
  train_data <- modelData[sample_indices, ]
  test_data <- modelData[-sample_indices, ]
  #Model
  ancestry_RFmodel <- randomForest(Ancestry ~ ., data = train_data, importance = TRUE, ntree = 100)
  #Accuracy
  predictions <- predict(ancestry_RFmodel, newdata = test_data)
  Accuracy <- mean(predictions == test_data$Ancestry)
  prediction <- data.frame(Predicted = predict(ancestry_RFmodel, newdata = samplePredVCF,  type = "prob"))  %>% gather(Ancestry, Probability, 1:3) %>% 
    mutate(Ancestry = gsub("Predicted.", "", Ancestry)) 
  if(max(prediction$Probability) < 0.60){
    prediction <- prediction %>%  arrange(desc(Probability)) %>% .[1:2,] %>% 
      summarise(Ancestry = paste(Ancestry, collapse = "/"), Probability = paste(Probability, collapse = "/"))
  }else{
    prediction <- prediction %>%  filter(Probability == max(Probability)) %>%
      summarise(Ancestry = paste(Ancestry, collapse = "/"), Probability = paste(Probability, collapse = "/"))
  }
  output <- prediction %>% 
    mutate(Info = VARIANTS, Accuracy = Accuracy, SampleID_Flowcell = VCF_list$SampleID_Flowcell[i]) %>% select(SampleID_Flowcell, everything())
  }else{
  output <- data.frame(SampleID_Flowcell = VCF_list$SampleID_Flowcell[i], Ancestry = "Insufficent data", Probability = 0, Info = 0, Accuracy = 0)
  }   
  print(output)
  fwrite(output, paste0(nxf_work_dir,"/../QC_",version,"/",VCF_list$Flowcell[i],"/",VCF_list$SampleID_Flowcell[i], ".ancestryPrediction.txt"), sep = "\t")
  
}
}

fwrite(VCF_list, "ancestry_samples.txt", sep = "\t")



# # #importance_df <- do.call("rbind", importance_list)
#  preds <- data.frame(pred = list.files(paste0(nxf_work_dir,"/QC/"), "*.ancestryPrediction.txt", recursive = TRUE, full.names = TRUE)) 
# # 
#  collect_preds <- list()
#  for(i in 1:nrow(preds)){
#    collect_preds[[i]] <- fread(preds$pred[i])
# #   
#  }
#  preds_df  <- do.call("rbind", collect_preds) 
#  
#  overview <- preds_df %>% mutate(Major = substr(Major, 0, 1)) %>% group_by(Major) %>% count(Ancestry)
