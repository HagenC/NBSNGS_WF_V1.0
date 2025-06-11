#!/usr/bin/env Rscript

#R-libraries:
#options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("data.table", "dplyr", "tidyr","stringr", "stringi", "magrittr", "DBI", "RSQLite", "randomForest")
for (pkg in packages) {
#  if (!require(pkg, character.only = TRUE)) {
#    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
#}

# for (pkg in packages) {
#   if (!require(pkg, character.only = TRUE)) {
#     install.packages(pkg,
#                      repos = "https://cloud.r-project.org",
#                      character.only = TRUE)
#     library(pkg, character.only = TRUE)
#   }
# }
# suppressPackageStartupMessages(suppressMessages(if (!require("data.table")) install.packages("data.table")))
# suppressPackageStartupMessages(suppressMessages(library("data.table")))
# suppressPackageStartupMessages(suppressMessages(if (!require("dplyr")) install.packages("dplyr")))
# suppressPackageStartupMessages(suppressMessages(library(dplyr)))
# suppressPackageStartupMessages(suppressMessages(if (!require("magrittr")) install.packages("magrittr")))
# suppressPackageStartupMessages(suppressMessages(library(magrittr)))
# suppressPackageStartupMessages(suppressMessages(if (!require("tidyr")) install.packages("tidyr")))
# suppressPackageStartupMessages(suppressMessages(library(tidyr)))
# suppressPackageStartupMessages(suppressMessages(if (!require("stringr")) install.packages("stringr")))
# suppressPackageStartupMessages(suppressMessages(library(stringr)))
# suppressPackageStartupMessages(suppressMessages(if (!require("stringi")) install.packages("stringi")))
# suppressPackageStartupMessages(suppressMessages(library(stringi)))
# suppressPackageStartupMessages(suppressMessages(if (!require("DBI")) install.packages("DBI")))
# suppressPackageStartupMessages(suppressMessages(library(DBI)))
# suppressPackageStartupMessages(suppressMessages(if (!require("RSQLite")) install.packages("RSQLite")))
# suppressPackageStartupMessages(suppressMessages(library(RSQLite)))
# suppressPackageStartupMessages(suppressMessages(if (!require("randomForest")) install.packages("randomForest")))
# suppressPackageStartupMessages(suppressMessages(library(randomForest)))

#Effect SnpEffect filter:
keep.SnpEff.EFFECTs <- c("conservative_inframe_deletion", "frameshift_variant", "splice_donor_variant", "stop_gained",
                         "disruptive_inframe_deletion", "missense_variant", "splice_acceptor_variant", "start_lost", "stop_lost",
                         "chromosome_number_variation", "exon_loss", "rare_amino_acid_variant", "transcript_ablation", "coding_sequence_variant",
                         "conservative_inframe_insertion", "disruptive_inframe_insertion", "regulatory_region_ablation", "TFBS_ablation")

#suppressPackageStartupMessages(suppressMessages(if (!require("rmarkdown")) install.packages("rmarkdown")))
#suppressPackageStartupMessages(suppressMessages(library(rmarkdown)))
#suppressPackageStartupMessages(suppressMessages(if (!require("doParallel")) install.packages("doParallel")))
#suppressPackageStartupMessages(suppressMessages(library(doParallel)))
#registerDoParallel(cores = 5)