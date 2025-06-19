# NBSNGS_WF_V1.0
NBSNGS_WF_V1.0 repository


1. Organization
NBSNGS_WF_V1.0/
├── README.md
├── main.nf                      # Main Nextflow workflow
├── nextflow.config             # Workflow parameters
├── condaEnvs/                  # Conda environments temp safe-house
├── envs/
│   └── environment.yaml        # Conda environment definitions
├── profiles/
│   └── profiles.config         # Configuration profiles
├── modules/                    # Module scripts (*.nf)
├── subworkflows/
│   ├── collectdata.nf
│   ├── mapping.nf
│   ├── calling.nf
│   └── reporting.nf
├── bin/                        # Executable pipeline scripts
│   └── <script>.R*
├── DATA/                       # Symbolic link to /srv/data/ILMN/RunFolder_NBSNGS/
├── assets/
│   ├── phenotype/
│   │   ├── FounderVariants.bed             # SSI focus variants
│   │   ├── geneFilter.list                 # Gene filter list for OPL
│   │   └── phenotype.gene.list
│   ├── sanity/
│   │   └── SnpEff_geneFilter.list_observation_count.tsv
│   ├── sql/
│   │   ├── SampleSheets.sqlite             # Sample sheet DB
│   │   └── variants.sqlite                 # Variant/QC DB
│   ├── temp/                               # Temporary files
│   └── tools/
│       └── vcfEffOnePerLine.pl
├── variants_/
│   ├── <sampleid_flowcell>.GATK.g.vcf.gz
│   ├── <sampleid_flowcell>.GATK.g.vcf.gz.tbi
│   ├── <sampleid_flowcell>.GATK.OPL.vcf
│   ├── <sampleid_flowcell>.lofreqDefault.vcf.gz
│   └── <sampleid_flowcell>.lofreqRaw.vcf.gz
├── QC_/
│   ├── multiqc_report.html
│   ├── <sampleid_flowcell>.ancestryPrediction.txt
│   ├── <sampleid_flowcell>.cov.gz
│   ├── <sampleid_flowcell>.mtDNAhg_classified.txt
│   ├── duplication_metrics/
│   │   └── <sampleid_flowcell>.txt
│   └── fastq/
│       └── <sampleid_flowcell>_fastqc.zip/.html
├── REPORT/
│   └── *.report.html
├── log_/                        # Workflow log summaries
│   ├── WF-reports.html
│   ├── haplotypecaller/
│   ├── lofreq/
│   └── mapping/

    ├── README.md
    ├── DATA  	                          #symbolic links to data folder (/srv/data/ILMN/RunFolder_NBSNGS/)     
    ├── assets                                #gets populated with  clinvarPathogenicTargetsubet.tsv, VOIS tables ets..     
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
    ├── bin                                   # Executable pipeline scripts
    │   └── <script>.R*
    ├── profiles
    │   └── profiles.config                   # Configuration profiles for compute environments
    ├── envs
    │   └── <name>/
    │       └── environment.yaml              # Conda environment definitions
    ├── main.nf                               # Main workflow 
    ├── modules/
    │   └── <module>.nf                       # Module scripts
    ├── subworkflows                          # Sub-workflows
    │   ├── collectdata.nf
    │   ├── mapping.nf
    │   ├── calling.nf
    │   └── reporting.nf
    ├── nextflow.config                       # Workflow parameters
    ├── condaEnvs                        #conda environments temporary safe-house.   
├───├
├─variants_<version>                        
│   └── <flowcell>                      
│            ├──<sampleid_flowcell>.GATK.g.vcf.gz/tbi          
│            ├──<sampleid_flowcell>.GATK.OPL.vcf  
│            ├──<sampleid_flowcell>.lofreqDefault.vcf.gz  
│            └──<sampleid_flowcell>.lofreqRaw.vcf.gz
├─ QC_<version>
│   └── <flowcell>                      
│           ├── *multiqc_report.html
│           └──<sampleid_flowcell>      
│                 ├──<sampleid_flowcell>.ancestryPrediction.txt
│                 ├──<sampleid_flowcell>.cov.gz
│                 ├──<sampleid_flowcell>.mtDNAhg_classified.txt
│                 ├──duplication_metrics
│                 │        └──<sampleid_flowcell>.txt
│                 └──fastq
│                      └──<sampleid_flowcell>_fastqc.zip/html
├─REPORT_<version>                       #Reports
│   └── <flowcell>
│          └── *.report.html
├─log_<version>                        #WF-reports.html and  mapping -, haplotypecaller & Lofreq logs.
     ├──<flowcell>
            ├──haplotypecaller
            ├──lofreq
            └──mapping
