# NBSNGS_WF_V1.0
NBSNGS_WF_V1.0 instructions:

## Clone this repository:
`git clone https://github.com/HagenC/NBSNGS_WF_V1.0.git .`

## How to setup nextflow:
  1. Install java in conda env:   
    `conda activate base`  
    `conda install openjdk=22.0.1`   (NGGT runs openjdk version "1.8.0_412"!)  
  2. Get nextflow (local):
     `$cd <path>/V1.0`  
     `wget -qO- https://get.nextflow.io | bash`  `   

## Quick start: 
  1.  `cd <path>/V1.0`  
  2.  `/nextflow run main.nf`


## Organization
<pre> <code> 
├─ README.md 
├───├── DATA                                    #symbolic links to data folder (/srv/data/ILMN/RunFolder_NBSNGS/)     
    ├── assets                                  #gets populated with  clinvarPathogenicTargetsubet.tsv, VOIS tables ets..     
    │   ├──phenotype
    │   │     ├──FounderVariants.bed            #SSI focus variant [CHR_POS][INFO][GENE][VARID]
    │   │     ├──geneFilter.list                #SnpEff target gene annotations list - used for filtering OPL 
    │   │     └──phenotype.gene.list     
    │   ├──sanity
    │   │     └──SnpEff_geneFilter.list_observation_count.tsv         #Target genes vs. SnpEff annoation check file
    │   ├──sql
    │   │     ├──SampleSheets.sqlite            #SampleSheet SQL database
    │   │     └──variants.sqlite                #Variant/QC SQL database
    │   ├──temp
    │   │     ├──temporary files      
    │   └──tools
    │        └──vcfEffOnePerLine.pl
    ├── bin                                     #Executable pipeline scripts
    │   └── < script >.R*
    ├── profiles
    │   └── profiles.config                     #Configuration profiles for compute environments
    ├── envs
    │   └── < name >
    │       └── environment.yaml                #Conda environment definitions
    ├── main.nf                                 #Main workflow 
    ├── modules/
    │   └── < module >.nf                       #Module scripts
    ├── subworkflows                            #Sub-workflows
    │   ├── collectdata.nf
    │   ├── mapping.nf
    │   ├── calling.nf
    │   └── reporting.nf
    ├── nextflow.config                         #Workflow parameters
    ├── condaEnvs                               #conda environments temporary safe-house.   
├───├
├─variants_< version >                        
│   └── < flowcell >                      
│            ├──< sampleid_flowcell >.GATK.g.vcf.gz/tbi          
│            ├──< sampleid_flowcell >.GATK.OPL.vcf  
│            ├──< sampleid_flowcell >.lofreqDefault.vcf.gz  
│            └──< sampleid_flowcell >.lofreqRaw.vcf.gz
├─ QC_< version >
│   └── < flowcell >                      
│           ├── *multiqc_report.html
│           └──< sampleid_flowcell >      
│                 ├──< sampleid_flowcell >.ancestryPrediction.txt
│                 ├──< sampleid_flowcell >.cov.gz
│                 ├──< sampleid_flowcell >.mtDNAhg_classified.txt
│                 ├──duplication_metrics
│                 │        └──< sampleid_flowcell >.txt
│                 └──fastq
│                      └──< sampleid_flowcell >_fastqc.zip/html
├─REPORT_< version >                            #Reports
│   └── < flowcell > 
│          └── .report.html
├─log_< version >                               #WF-reports.html and  mapping -, haplotypecaller & Lofreq logs.
     ├──< flowcell >
            ├──haplotypecaller
            ├──lofreq
            └──mapping
        </code> </pre>

---

## Crtitical process steps

#### Process: SampleSheet Collection
- Searches for `SampleSheet.csv` i `$WD/DATA/*`
- Verifies presence of accompanying `CopyComplete.txt`  `$WD/DATA/*`
- Imports only `SampleSheet.csv` files containing all required columns:  
  `["Sample_ID", "Sample_Name", "Sample_Project", "Description"]`
- Validates Flowcell ID using regex pattern:  
  `[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+`
- Queries only samples with `Project_ID` matching any of:  
  `"NBS-NGS"`, `"NBS_NGS"`, `"nbs-ngs"`

#### Process: Build Variant SQLdatabase
- Only gene annotations listed in [Gene filter list](v1.0/assets/phenotype/geneFilter.list) are considered for reporting

#### Process: Report Generation
- Description/phenotype selection must have a complete match with `Phenotype` in  [Phenotype - Gene](v1.0/assets/phenotype/phenotype.gene.list)

---

## Support data information

- **Reference Genome**:  
  - [Homo_sapiens_assembly38.fasta](https://storage.cloud.google.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta)  
  - Masked as described in: [DOI: 10.1038/s41587-021-01158-1](https://doi.org/10.1038/s41587-021-01158-1)  
  - Alternative scaffolds removed

- **GATK Base Quality Score Recalibration (BQSR) Files**:
  - [Mills + 1000G Gold Standard Indels (hg38)](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz)
  - [1000G Phase 1 High Confidence SNPs (hg38)](https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz)


- **SnpEff Database**: `GRCh38.113` (custom/homebrew build)
  - **Gene annotation (GTF)**:  
    (ftp://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz)
  - **Reference genome (FASTA)**:  
    (ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)

- **gnomAD**: `Genomes v.4.1`  
    https://gnomad.broadinstitute.org/data

 - **Ancestry reference data**:  
  [1000genomes](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/)

- **xGen™ Exome Hybridization Panel**:
  [target bed file](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/supplementary-product-info/xgen-exome-hyb-panel-v2-targets-manifest-hg38.txt?sfvrsn=b6c0e207_2)






