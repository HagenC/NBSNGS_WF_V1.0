## Detailed description:
### Workflow Overview: `main.nf`

**Subworkflows:**

- [`COLLECTDATA`](#collectdata)
- [`MAPPING`](#mapping)
- [`CALLING`](#calling)
- [`ANNOTATION`](#annotation)
- [`REPORTING`](#reporting)

### COLLECTDATA [`collectdata.nf`]

**Processes:** - `FINDDATA [finddata.nf]`  
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
  - **Input:** `DOWNSAMPLE.out` (dummy)  
  - **Tool:** `fastqc=0.12.1`

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

---
## Detailed R scripts:
---

### COLLECTDATA [`collectdata.nf`]

**Subworkflow Process:** `FINDDATA` [`finddata.nf`]  
- **Tool:** `R 4.3.1`  
- **Script:** `sqlsamplesheet.R`


#### 1. SampleSheet Collection
- Searches for `SampleSheet.csv` in `$WD` (**CRITICAL**)
- Verifies presence of accompanying `CopyComplete.txt` in `$WD` (**CRITICAL**)
- Collects only those `SampleSheet.csv` files that have a matching `CopyComplete.txt` (**CRITICAL**)


#### 2. SampleSheet Import and Validation
- Imports only `SampleSheet.csv` files containing all required columns:  
  `["Sample_ID", "Sample_Name", "Sample_Project", "Description"]`
- Validates Flowcell ID using regex pattern:  
  `[[:digit:]]{6}_[[:alpha:]]+[[:digit:]]+_[[:digit:]]+_[[:alnum:]]+` (**CRITICAL**)
- Generates unique `SampleID.Flowcell` identifier: `Sample_ID + Flowcell`
- Removes entries with missing `SampleID`, missing `Flowcell`, or duplicate rows


#### 3. FASTQ File Discovery
- Recursively searches collected paths for `*.fastq.gz` files
- Constructs `SampleID.Flowcell` identifiers and filters for matches in the imported SampleSheet (**CRITICAL**)
- Captures `filesize` and `mtime` for each valid FASTQ file


#### 4. Merging FASTQ and SampleSheet Data
- Combines FASTQ metadata with SampleSheet entries
- Removes SampleSheet entries that lack corresponding FASTQ files


#### 5. SQLite Database Construction
- Database path: `@WD/assets/SampleSheets.sqlite`
- Checks for existing table `SAMPLESHEETS`
  - If absent, creates the table
  - If present, appends new entries where `SampleID.Flowcell` is not already listed


#### 6. Sample Selection for Processing
- Queries samples with `Project_ID` matching any of:  
  `"NBS-NGS"`, `"NBS_NGS"`, `"nbs-ngs"` (**CRITICAL**)

#### 7. Processed Sample Check
- Skips samples where `SampleID_Flowcell.GATK.g.vcf.gz` already exists

#### 8. Processing Prioritization
- Prioritizes samples based on latest `mtime`
- Excludes already processed samples
- Limits selection to **11 samples max**


#### 9. Outputs
- `Sample.collection.info`: Existing samples from `SampleSheets.sqlite` with relevant `Project_ID`s
- `ProcessedSamples.txt`: Sent via email upon successful completion
- `ProcessingIDs.txt`: Input to downstream processes (`channel emit: alignment_tuple`)

		

### Process: `MTDNA_HAPLOGROUP` [`mtdna_haplogroup.nf`]

- **Input:** `OPL_ch`  
- **R Script:** `mtdnahaplogroup.R`
- **Tools:** `R=4.3.1`, `haplogrep=2.4.0`

---

#### 1. Search Existing Results
- Searches for existing classification files:  
  `*.mtDNAhg_classified.txt` in `/QC/`

#### 2. Filter New Samples
- Looks for `*.g.vcf.gz` files in `/variants/`
- Filters out `SampleID_Flowcell` entries that already have corresponding classification files


#### 3. Predict mtDNA Haplogroup
    - classifies `g.vcf.gz`using haplogrep.
    -  **Condition:** Classifications successfull?
    - **If TRUE:**
    - Computes QC metrics: mean depth (DP), variant count
    - **If FALSE:**
    - Returns default placeholder.
---

#### 4. Output
- Writes classification files to:  
  `/QC/*.mtDNAhg_classified.txt`

---

### Process: `SEX` [`sex.nf`]

- **Input:** `OPL_ch`, `bed`  
- **R Script:** `sexcheck.R`
- **Tools:** `R=4.3.1`

#### 1. Find Input Files
- Locates:
  - Coverage files: `/QC/*.cov.gz`
  - Existing sex prediction results: `/QC/*.predictedSex.txt`
- Filters out already classified `SampleID_Flowcell`

#### 2. Predicting Sex
- **Condition:** Number of `chrX` coverage entries > 1
  - **If TRUE:**
    - Computes mean coverage of `chrX` and `chrY`
    - **Rule:**  
      - If `mean(chrX)/mean(chrY) < 3`, predict `"Male"`  
      - Else, predict `"Female"`
  - **If FALSE:**
    -- Returns default placeholder. 


### Process: `ANCESTRY` [`ancestrymodel.nf`]

- **Input:** `OPL_ch`  
- **R Script:** `ancestrymodel.R`  
- **Channel Emit:** `ANCESTRY.out`

#### 1. Extract Sample Metadata
- Loads `SAMPLESHEETS` from SQLite:  
  `/assets/SampleSheets.sqlite`
- Filters for `Sample_Project` in:  
  `"NBS-NGS"`, `"nbs-ngs"`, `"NBS_NGS"`

#### 2. Find Input Files
- Existing results: `/QC/*.ancestryPrediction.txt`
- New candidate VCFs: `/variants/*GATK.vcf.gz`  
  - Filters out `SampleID_Flowcell` values already classified

#### 3. Ancestry Prediction

##### Reference Data Imports:
1. Selected RS IDs:  
   `parmas.ancestry_RS_selection`
2. Reference metadata:  
    `ancestry_reference ` 

##### Sample Data:
4. Imports target sample VCF  
5. Computes QC measure:  
   - Count of variants in sample VCF that match reference RS IDs

---

#### 4. Classification Logic

- **Condition:** `variants > 100`
  - **If TRUE:**
    1. Genotypes (GT) are matched and converted to integers
    2. Data is reshaped to long format
    3. Reference metadata is subset to matching RS IDs
    4. Splits data: 80% train / 20% test
    5. Trains `randomForest` model, tests on held-out data
    6. Predicts ancestry for the sample
    7. Reports:
       - If prediction probability < 0.6 → **Top 2 predictions**
       - Else → **Top 1 prediction**
       
  - **If FALSE:**
    - Returns default placeholder.

---

Let me know if you'd like a version with tables, diagrams, or callouts for key thresholds like `INFO > 100` or `probability < 0.6`.


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
