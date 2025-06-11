#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
nxf_work_dir <- args[1]
version <- args[2]
#nxf_work_dir <-  "/srv/data/ILMN/PIPELINES/NBSNGS_WF_v1.0/v1.0/"
#version <- "v.1.0"
source(paste0(nxf_work_dir,"/bin/dependencies.R"))


existing_mtDNAhg <- data.frame(mtDNA = list.files(paste0(nxf_work_dir,"/../QC_",version,"/"),"*.mtDNAhg_classified.txt", recursive = TRUE, full.names = TRUE)) %>% 
  mutate(SampleID_Flowcell = str_sub(mtDNA, stri_locate_last_fixed(mtDNA, "/")[,1]+1, nchar(mtDNA))) %>%
  mutate(SampleID_Flowcell  = sub(".mtDNAhg_classified.txt", "", SampleID_Flowcell ))  
VCFs <- data.frame(VCF = list.files(paste0(nxf_work_dir,"/../variants_",version,"/"), "*.GATK.g.vcf.gz$", recursive = TRUE, full.names = TRUE)) %>%
  mutate(SampleID_Flowcell = str_sub(VCF,stri_locate_last_fixed(VCF, "/")[,1]+1, nchar(VCF))) %>%
  mutate(SampleID_Flowcell  = sub(".GATK.g.vcf.gz", "", SampleID_Flowcell )) %>% 
  mutate(Flowcell = str_extract(SampleID_Flowcell, "[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+")) %>% filter(!SampleID_Flowcell %in% existing_mtDNAhg$SampleID_Flowcell ) 

i <- 2
if(nrow(VCFs) >= 1){
for(i in 1:nrow(VCFs)){ 
 system(paste0("haplogrep classify  --format=vcf --input=",VCFs$VCF[i], " --output=",nxf_work_dir,"/../QC_",version,"/",VCFs$Flowcell[i],"/",VCFs$SampleID_Flowcell[i],".mtDNAhg_classified.txt --extend-report"))
 mtDNA_prediction <- fread(paste0(nxf_work_dir,"/../QC_",version,"/",VCFs$Flowcell[i],"/",VCFs$SampleID_Flowcell[i],".mtDNAhg_classified.txt")) %>% select(-Remaining_Polys)
 print(nrow(mtDNA_prediction))
 if(nrow(mtDNA_prediction) == 1){
     mtDNA_QC <- fread(cmd = paste0("zcat ", VCFs$VCF[i], " | grep -v ^##")) %>% filter(`#CHROM` == "chrM") %>% separate(10, c("GT","AD", "DP"), sep = ":") %>% 
       filter(GT != "0/0") %>% 
       mutate(DP = as.numeric(DP)) %>%  summarise(Median.DP = median(DP, na.rm = TRUE), Variants = n())
     mtDNAhg_classified.txt <- cbind(mtDNA_prediction, mtDNA_QC)
  }else{
   mtDNAhg_classified.txt <- data.frame(SampleID = VCFs$SampleID_Flowcell[i]) %>%   
                             mutate( Haplogroup = "Insufficient data",
                                     Rank = "0",
                                     Quality = "0",
                                     Range = "0",
                                     Not_Found_Polys = "-",
                                     Found_Polys = "-",
                                     AAC_In_Remainings = "-",
                                     Input_Sample = "-",
                                     Median.DP = "0",
                                     Variants = "0")
   
    }
 fwrite(mtDNAhg_classified.txt,paste0(nxf_work_dir,"/../QC_",version,"/",VCFs$Flowcell[i],"/",VCFs$SampleID_Flowcell[i],".mtDNAhg_classified.txt"), sep = "\t")
  }
}  

fwrite(existing_mtDNAhg, "mtDNA_haplotypedSamples.txt", sep = "\t")
