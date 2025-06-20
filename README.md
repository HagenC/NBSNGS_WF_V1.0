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


## Detailed description:
### Workflow Overview: `main.nf`

**Subworkflows:**

- [`COLLECTDATA`](#collectdata)
- [`MAPPING`](#mapping)
- [`CALLING`](#calling)
- [`ANNOTATION`](#annotation)
- [`REPORTING`](#reporting)

### COLLECTDATA [`collectdata.nf`]

**Processes:**

- `FINDDATA [finddata.nf]`  
  - **Tool:** `R=4.3.1`  
  - **Script:** `sqlsamplesheet.R`  
  - **Function:**  
    - Finds and validates `SampleSheet.csv` with `CopyComplete.txt`  
    - Combines metadata with FASTQ files  
    - Creates SQLite DB: `SampleSheets.sqlite`  
    - Emits channel: `alignment_tuple`  

---

### MAPPING [`mapping.nf`]

**Input Channels:** `samples_for_alignment`, `unique_flowcells` created from `alignment_tuple`

**Processes:** `DOWNSAMPLE [downsample.nf]`  
  - **Tool:** `bbmap=39.08`  
  - **Settings:** `samplereadstarget=93487190`  

- **Process:** `FASTQC [fastqc.nf]`  
  - **Input:** `DOWNSAMPLE.out`  
  - **Tool:** `multiqc=1.22.2`

- **Channel:** `FASTQC.out`

- **Process:** `ALIGNMENT [alignment.nf]`  
  - **Input:** `DOWNSAMPLE.out`  
  - **Tools:** `bwa-mem2=2.2.1`, `samtools=1.20`  
  - **External Input:** `params.bed`

- **Process:** `RAW_INDEX [index.nf]`  
  - **Input:** `ALIGNMENT.out`  
  - **Tool:** `samtools=1.20`

- **Process:** `ADDREADGROUP [addreadgroup.nf]`  
  - **Input:** `RAW_INDEX.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `MARKDUPLICATES [markduplicates]`  
  - **Input:** `ADDREADGROUP.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `MULTIQC [multiqc.nf]`  
  - **Input:** `FASTQC.out`, `MARKDUPLICATES.out`  
  - **Tool:** `multiqc=1.22.2`  
  - **Conditional Execution:** `doQC = FALSE` (default)

<details>
<summary>🔹 If <code>doQC = TRUE</code></summary>

- **Process:** `RAW_DEPTH [mosdepth.nf]`  
  - **Input:** `RAW_INDEX.out`  
  - **Tool:** `mosdepth=0.3.3`

- **Process:** `HSMETRICS [hs_metrics.nf]`  
  - **Input:** `ALIGNMENT.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `ALIGNMENT_METRICS [alignment_metrics.nf]`  
  - **Input:** `ALIGNMENT.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `INSERT_SIZE_METRICS [insert_size_metrics.nf]`  
  - **Input:** `ALIGNMENT.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Channel:** `FASTQC.out`, `HSMETRICS.out`, `ALIGNMENT_METRICS.out`, `INSERT_SIZE_METRICS.out`, `RAW_DEPTH.out`

- **Process:** `MULTIQC [multiqc.nf]`  
  - **Input:** `FASTQC.out`, `HSMETRICS.out`, `ALIGNMENT_METRICS.out`, `INSERT_SIZE_METRICS.out`, `RAW_DEPTH.out`  
  - **Tool:** `multiqc=1.22.2`

</details>

- **Channel Emit:** `INDEX.out`

### CALLING [`calling.nf`]
- **Conditional Execution:** `doCalling = TRUE` (default)

- **Process:** `CLEANUP [cleanup.nf]`  
  - **Input:** `MAPPING.emit`  
  - **Tool:** `samtools=1.20`

- **Process:** `INDELQUAL [indelqual.nf]`  
  - **Input:** `CLEANUP.out`  
  - **Tool:** `lofreq=2.1.5`

- **Process:** `INDEX_CLEAN [index.nf]`  
  - **Input:** `INDELQUAL.out`  
  - **Tool:** `lofreq=2.1.5`

- **Process:** `LOFREQNODEFAULT [lofreq_nodefault.nf]`  
  - **Input:** `INDEX_CLEAN.out`  
  - **Tool:** `lofreq=2.1.5`

- **Process:** `LOFREQCALL [lofreq_default.nf]`  
  - **Input:** `INDEX_CLEAN.out`  
  - **Tool:** `lofreq=2.1.5`

- **Process:** `BQSR [bqsr.nf]`  
  - **Input:** `CLEANUP.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `APPLY_BQSR [apply_bqsr.nf]`  
  - **Input:** `BQSR.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `INDEX_HAPLOTYPECALLER [index.nf]`  
  - **Input:** `APPLY_BQSR.out`  
  - **Tool:** `samtools=1.20`

- **Process:** `COVERAGE [coverage.nf]`  
  - **Input:** `INDEX_HAPLOTYPECALLER.out`  
  - **Tool:** `samtools=1.20`

- **Process:** `HAPLOTYPECALLER [haplotypecaller.nf]`  
  - **Input:** `INDEX_HAPLOTYPECALLER.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `GTOVCF [gtovcf.nf]`  
  - **Input:** `HAPLOTYPECALLER.out`  
  - **Tool:** `gatk4=4.5.0.0`

- **Channel Emit:** `GTOVCF.out`

### ANNOTATION [`annotating.nf`]
- **Conditional Execution:** `doReporting = TRUE` (default)

- **Process:** `UPDATE_CLIVAR [update_clinvar.nf]`  
  - **Input:**  
    - `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz/tbi`  
    - `params.bed`  
    - `params.OPL`  
  - **Output:** `/assets/clinvarPathogenicTargetsubet.tsv`  
  - **Tools:** `bash=4.2.46`, `bcftools=1.6`, `SnpSift`

- **Process:** `GATK_LA [leftAlign.nf]`  
  - **Input:** `CALLING.emit`  
  - **Tool:** `gatk4=4.5.0.0`

- **Process:** `ANNOTATE_CLINVAR [annotate_clinvar.nf]`  
  - **Input:** `GATK_LA.out`  
  - **Tools:** `snpsift=5.2/snpeff=5.2`, `clinvar.vcf.gz`

- **Process:** `ANNOTATE_GNOMAD [annotate_gnomad.nf]`  
  - **Input:** `ANNOTATE_CLINVAR.out`  
  - **Tools:** `snpsift=5.2/snpeff=5.2`, `gnomad.4.1`

- **Process:** `ANNOTATE_ONTOLOGY [annotate_ontology.nf]`  
  - **Input:** `ANNOTATE_GNOMAD.out`  
  - **Tools:** `snpsift=5.2/snpeff=5.2`

- **Process:** `OPL [opl.nf]`  
  - **Input:** `ANNOTATE_ONTOLOGY.out`  
  - **Tools:** `snpsift=5.2/snpeff=5.2`

- **Channel:** `OPL_ch = OPL.out`

- **Process:** `MTDNA_HAPLOGROUP [mtdna_haplogroup.nf]`  
  - **Input:** `OPL_ch`  
  - **Tools:** `haplogrep3`, `R=4.3.1`

- **Process:** `SEX [sex.nf]`  
  - **Input:** `OPL_ch`  
  - **Tool:** `R=4.3.1`

- **Process:** `ANCESTRY [ancestrymodel.nf]`  
  - **Input:** `OPL_ch`  
  - **Tool:** `R=4.3.1`

- **Channel Emit:** `ANCESTRY.out`

 ### REPORTING [`reporting.nf`]
- **Conditional Execution:** `doReporting = TRUE` (default) 
 **Input dummy:**  `ANNOTATION.emit`

- **Process:** `BUILDSQL [buildsql.nf]`  
  - **Tool:** `R=4.3.1`
  - **Output:** `variants.sqlite` 

- **Process:** `VOIDATABASE [voidatabase.nf]`  
  - **Tool:** `R=4.3.1`
  - **Output:** `VOI_DB`  
  - **Output:** `VOI_DB_PHENOTYPE`  

- **Process:** `RENDER [render.nf]`  
  - **Tool:** `R=4.3.1`
  - **Output:** `*< sampleid_flowcell_description > .report.html`  



### Process: `MTDNA_HAPLOGROUP` [`mtdna_haplogroup.nf`]
**Input:** `OPL_ch`  
**R Script:** `mtdnahaplogroup.R`

1. **Search existing results**
   - Find all `*.mtDNAhg_classified.txt` files in `/QC/`

2. **Filter new samples**
   - Search for `*.GATK.OPL.vcf` in `/variants/`
   - Exclude `SampleID_Flowcell` values that already have a classification file

3. **Predict mtDNA haplogroup**
- process: MTDNA_HAPLOGROUP [mtdna_haplogroup.nf]; input:  OPL_ch
	-mtdnahaplogroup.R:
		1. Finding existing "/QC/" *.mtDNAhg_classified.txt"
		2. Finding OPLs in "/variants" - "*.GATK.OPL.vcf$" and filtering out existing SampleID_Flowcell from step 1. 
		3. Predicting mtDNA haplogroup:
			1. Conditional: Checking if "chrM" exists in OPL 
				1. TRUE: Importing OPL CHROM = "chrM" and extracting GT == "1/1"
				2.Creating "QC" - mean DP and variant count. 
				3. Conditional: Checking if number of variants in step 1 is > 1.
					1.Filter for SNPs and creating posAlt variable.
					2. Writing out step 1.
					3. Using haplogrep3 to predict "classify --tree=phylotree-rcrs@17.2"
					4. Adding "QC" to prediction and changign predicted Haplogroup to "Insufficient coverage" if mean DP is < 20. 
				2.FALSE: Creating with data.frame(SampleID = OPL$SampleID_Flowcell[i], Haplogroup = "Insufficient data", Rank = 0, Quality = 0, Range = 0-0, Major = "-", 
                                    PredictedAncestry = "-",  Mean.DP = 0,  n = 0)
		4. OUTPUT: "/QC/" *.mtDNAhg_classified.txt"
- process: SEX [sex.nf]; input: OPL_ch, bed, bed-padded.
	-sexchech.R:
		1. Conditional: if(file.exists(paste0(nxf_work_dir,"/assets/bedFile_expanded.bed"))== "FALSE")
			FALSE: Import bed-file and expand to chr:pos to make expanded-bed-file.
			TRUE: import expanded-bed-file.
		2.FIND files:
			1. COVERAGE "/QC/"  "*.cov.gz$
			2. SEX "/QC/" *.predictedSex.txt" filtering out SampleID_Flowcell in step 1. 
		2. Predicting sex:
			1. Conditional number of lines in COVERAGE CHROM == chrX
				1. >1 mean of chrX and chrY  if ratio of mean chrX/chrY < 3 = male else female. 
				2. < 1 creating data.frame data.frame(chrX = 0, chrY = 0)%>% mutate( SampleID_Flowcell = DB$SampleID_Flowcell[i]) %>% 
					mutate(XY.ratio = chrX/chrY) %>% mutate(Sex.prediction = ifelse(XY.ratio < 3, "Male", "Female")) %>% 
					mutate(Sex.prediction = ifelse(is.na(Sex.prediction), "Insufficent data", Sex.prediction))
- process: ANCESTRY [ancestrymodel.nf]; input: OPL_ch
	- ancestrymodel.R 
		1.Extract SAMPLESHEETS from "/assets/SampleSheets.sqlite"
		2.Sample_Project = NBS-NGS, nbs-ngs or NBS_NGS 
		3. FIND files:
			DONE "QC/" "*.ancestryPrediction.txt"
			VCF_list "/variants" "*GATK.vcf.gz$" filtering out SampleID_Flowcell in step 1. 
		4. Ancestry prediction:
			1. Importing selected-RS from "/assets/FeatureSelection_MDA4.txt" (how is it made?)
			2. Importing reference-VCF "/assets/ancestry_reference_population.vcf" (how is it made?)
			3. Importing reference-DATA "/assets/ReferenceData.tsv" (how is it made?)
			4. Importing Sample-VCF
			5. Creating INFO as number of variants in in Sample-VCF that have matching variant in reference-VCF
			6. Conditional: INFO > 100
				1. If TRUE: GTs in Sample-VCF is matched with Reference-VCF IDs and GTs converted to integer and wrangled into long format. 
				2. Reference-DATA is subet to matching Sample-VCF RS ids and split into train and test set (0.8)
				3. randomforest model is build and tested on test set. 
				4. Ancestry is predicted on Sample-VCF subset adn Accuracy of test-set calculated.
				5. If Probability is < 0.6 top 2 is shown else top 1. 
			7. If FALSE: data.frame is generated: data.frame(SampleID_Flowcell = VCF_list$SampleID_Flowcell[i], Ancestry = "Insufficent data", Probability = 0, Info = 0, Accuracy = 0)
- channel emit: ANCESTRY.out




- conditional invocation: doReporting = TRUE default
Subworkflow: REPORTING [reporting.nf] ; input: ANNOTATION.emit, rmd-file [report.Rmd]
- process: BUILDSQL [buildsql.nf] ; input : bedfile, bedfile.padded ; env/tool: R=4.3.1, output: "SAMPLESforREPORTING.txt"
	- buildsql.R:
		1.Extract SAMPLESHEETS from "/assets/SampleSheets.sqlite"
		2.Sample_Project = NBS-NGS, nbs-ngs or NBS_NGS 
		3. Conditional: if(file.exists(paste0(nxf_work_dir,"/assets/bedFile_expanded.bed"))== "FALSE")
			FALSE: Import bed-file and expand to chr:pos to make expanded-bed-file.
			TRUE: import expanded-bed-file.
		4.Import "phenotype.gene.list"
		5.Filter expanded-bed-file for target gene-ids.
		6. FIND files:
			1. OPLs in "/variants" - "*.GATK.OPL.vcf$"
			2. LOFREQ in "/variants" - "*.lofreqDefault.vcf.gz$"
			3. COV in "/QC" - "*.cov.gz$"
			4. RAW in "/variants" - "*.lofreqRaw.vcf.gz$"
			5. SEX in "/QC" - *.predictedSex.txt$"
			6: mtDNA in "/QC/" - "*.mtDNAhg_classified.txt$"
			7: ANCESTRY  in "QC" - "*.ancestryPrediction.txt$"
		7. Creating sql-database in @WD/assets/variants.sqlite:
			1.Check if tables: LOFREQ, OPL, COV, RAW, SEX, mtDNA and ANCESTRY exists.
			2. TRUE table exists:
				1.For each table LOFREQ, OPL, COV, RAW, SEX, mtDNA and ANCESTRY extract SampleID_Flowcell and filter out SampleID_Flowcell from step 6.(FIND files:)
				2.Add if-any to the table.
			3. FALSE table exist:
				1. Create table and add to the table.
		8. OUTPUT: step 1. NBSNGS_Samplesheet where SampleID_Flowcell *".report.html" is filtered out. 
- process: VOIDATABASE [voidatabase.nf] ; input: BUILDSQL.out ; env/tool: R=4.3.1, output: VOI_DB, VOID_DB_PHENOTYPE
	- voiddatabase.R:
		1.Extract SAMPLESHEETS from "/assets/SampleSheets.sqlite"
		2.Sample_Project = NBS-NGS, nbs-ngs or NBS_NGS
		3.Import "phenotype.gene.list"
		4.Import FounderVariants.bed
		5.Importing SnpEff "EFFECT" filter:
		6.Importing from  "/assets/variants.sqlite" OPL table where EFFECT in EFFECT-filter and CLNSIG "*pathogen*" and not "*benign*"  - filtered for AF_ < 0.01 and "*Pathogenic"
		7.Creating chr:pos extraction list based on step 6. 
		8.Importing coverage from "/assets/variants.sqlite" COV table where CHR_POS in step 7.
		9.Creating "VOI_DB":
			1. Creating "coverage_VOIS_count" number of failed.alleles (DP < 10) and passed.alleles (PD >= 10)
			2.Creating "VOI_DB_GENE" VARID (chr-pos-ref-alt) and GENE from OPL_VOIS from step 6.
			3.Creating "VOI_DB" - VARID (chr-pos-ref-alt) GT (0/0,0/1,1/1) count caclulating AF_NBSNGS based on step 1 and 2. 
		10.Creating "VOI_DB_PHENO":
			1. "coverage_PHENO" number of failed.alleles (DP < 10) and passed.alleles (PD >= 10) geouped by "Description"
			2. Creating "VOI_DB_PHENOTYPE" - VARID (chr-pos-ref-alt), Description ,GT (0/0,0/1,1/1) count caclulating AF_NBSNGS based on step 1 grouped by Description.
- process: RENDER [render.nf]; input: VOIDATABASE.out ; env/tool: R=4.3.1
