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

- **Processes:** `DOWNSAMPLE [downsample.nf]`  
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

---

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

---

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

---

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
- Searches for `SampleSheet.csv` i `$WD/DATA/` (**CRITICAL**)
- Verifies presence of accompanying `CopyComplete.txt`  `$WD/DATA/` (**CRITICAL**)
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
- Database path: `@WD/assets/sql/SampleSheets.sqlite`
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


#### 1. Search Existing Results
- Searches for existing classification files:  
  `*.mtDNAhg_classified.txt` in `$WD/../QC/`

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

#### 4. Output
- Writes classification files to:  
  `$WD/../QC/*.mtDNAhg_classified.txt`

---

### Process: `SEX` [`sex.nf`]

- **Input:** `OPL_ch`, `bed`  
- **R Script:** `sexcheck.R`
- **Tools:** `R=4.3.1`

#### 1. Find Input Files
- Locates:
  - Coverage files: `$WD/../QC/*.cov.gz`
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
  `$WD/assets/sql/SampleSheets.sqlite`
- Filters for `Sample_Project` in:  
  `"NBS-NGS"`, `"nbs-ngs"`, `"NBS_NGS"`

#### 2. Find Input Files
- Existing results: `$WD/../QC/*.ancestryPrediction.txt`
- New candidate VCFs: `$WD/../variants/*GATK.vcf.gz`  
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

### Process: `BUILDSQL` [`buildsql.nf`]
 
**Tool:** `R 4.3.1`  
**Output:** `SAMPLESforREPORTING.txt`, `variants.sqlite` 
**R Script:** `buildsql.R`


#### 1. Load Sample Metadata
- Extracts `SAMPLESHEETS` from `$WD/assets/sql/SampleSheets.sqlite`
- Filters for `Sample_Project` in: `"NBS-NGS"`, `"nbs-ngs"`, `"NBS_NGS"`


#### 2. Load Gene Targets
- Imports: `phenotype.gene.list`
- Filters expanded BED file using gene IDs from phenotype list

#### 4. Locate Input Files
- Searches for:
  1. OPL VCFs: `$WD/../variants/*.GATK.OPL.vcf`
  2. LOFREQ VCFs: `$WD/../variants/*.lofreqDefault.vcf.gz`
  3. Coverage files: `$WD/../QC/*.cov.gz`
  4. Raw LOFREQ VCFs: `$WD/../variants/*.lofreqRaw.vcf.gz`
  5. Sex prediction results: `$WD/../QC/*.predictedSex.txt`
  6. mtDNA haplogroup classifications: `$WD/../QC/*.mtDNAhg_classified.txt`
  7. Ancestry predictions: `$WD/../QC/*.ancestryPrediction.txt`

#### 5. Build SQLite Database
- Creates or updates: ``$WD/assets/sql/variants.sqlite`

##### Table Handling:
- **Checks for existence of tables:**  
  `LOFREQ`, `OPL`, `COV`, `RAW`, `SEX`, `mtDNA`, `ANCESTRY`
- **If table exists:**
  - Extracts `SampleID_Flowcell` from input files
  - Filters out existing entries
  - Appends new data to the table
- **If table does not exist:**
  - Creates the table
  - Adds corresponding data

#### 6. Output
- Generates a filtered `NBSNGS_Samplesheet`
- Final output: `SAMPLESforREPORTING.txt`  
  - Contains samples where `*.report.html` does **not** already exist for `SampleID_Flowcell`

---

### Process: `VOIDATABASE` [`voidatabase.nf`]

**Input:** `BUILDSQL.out`  
**Tool/Environment:** `R 4.3.1`  
**Output:** `VOI_DB`, `VOID_DB_PHENOTYPE`  
**R Script:** `voiddatabase.R`

---

#### 1. Load Metadata
- Extracts `SAMPLESHEETS` from `@WD/assets/sql/SampleSheets.sqlite`
- Filters for `Sample_Project` in: `"NBS-NGS"`, `"nbs-ngs"`, `"NBS_NGS"`

#### 2. Load Reference Data
- Imports:
  - `phenotype.gene.list`
  - `FounderVariants.bed`
  - SnpEff `EFFECT` filter list

#### 3. Import and Filter Variants
- From `$WD/assets/sql/variants.sqlite`, `OPL` table:
  - Filters where:
    - `EFFECT` is in EFFECT filter list
    - `CLNSIG` contains `"*pathogen*"` and **not** `"*benign*"`
    - `AF_` < 0.01
    - Text contains `"Pathogenic"`

#### 4. Generate Variant Coordinates
- Creates `chr:pos` extraction list based on filtered variants from step 3


#### 5. Import Coverage Data
- From `$WD/assets/sql/variants.sqlite`, `COV` table:
  - Filters where `CHR_POS` matches the list from step 4

#### 6. Create `VOI_DB`
- **coverage_VOIS_count**:
  - Counts of failed alleles (`DP < 10`) and passed alleles (`DP ≥ 10`)
- **VOI_DB_GENE**:
  - Contains `VARID` (`chr-pos-ref-alt`) and `GENE` for variants from step 3
- **VOI_DB**:
  - Contains:
    - `VARID` (`chr-pos-ref-alt`)
    - Genotype counts (`0/0`, `0/1`, `1/1`)
    - Calculated `AF_NBSNGS` from the filtered sample data (steps 1 & 2)

#### 7. Create `VOI_DB_PHENOTYPE`
- **coverage_PHENO**:
  - Counts of failed and passed alleles grouped by `Description`
- **VOI_DB_PHENOTYPE**:
  - Contains:
    - `VARID` (`chr-pos-ref-alt`)
    - `Description`
    - Genotype counts (`0/0`, `0/1`, `1/1`)
    - Calculated `AF_NBSNGS` grouped by `Description`

--- 

### Process: `RENDER` [`render.nf`]

**Input:** `VOIDATABSE.out`  
**Tool/Environment:** `R 4.3.1`  
**Output:** `Individual HTML reports per sample (`*.report.html`)`, `Processed_SampleIDs.txt`  
**R Script:** `render.nf`

#### 1. Load Metadata
- Connects to:
  - `$WD/../assets/sql/SampleSheets.sqlite` → fetches samples with `Sample_Project` in `["NBS-NGS", "nbs-ngs", "NBS_NGS"]`
  - `variants.sqlite` → retrieves distinct `SampleID_Flowcell` from the `OPL` table
- Determines completed reports by scanning `REPORTS_<version>/` for `*.report.html`

#### 3. Identify Samples for Reporting
- Filters for samples with:
  - Entries in the `OPL` table
  - **No** existing report generated

#### 4. Import temp-files
- ClinVar pathogenic subset:  
  `clinvarPathogenicTargetsubet.tsv`
- BED files:  
  `bedfile_target.bed`, `bedfile_coverage_expanded.bed`, `bedfile_call_expanded.bed`
- Phenotype gene list:  
  `phenotype.gene.list`
- Optional founder variants:  
  `FounderVariants.bed`
- VOI population data:  
  `VOI_DB` and `VOI_DB_PHENOTYPE`

#### 5. Per Sample Processing (loop)
##### 5.1 Sample Info
- `SampleID_Flowcell`, `Description`, `Flowcell` (phenotype ID)

##### 5.2 Target Genes
- Filters `phenotype.gene.list` using the phenotype description

##### 5.3 Variant Annotation
- Queries `OPL`, `RAW`, `LOFREQ` tables for the sample
- Filters using:
  - SnpEff effect matches
  - CLNSIG (must be pathogenic, not benign)
  - `AF_nfe < 0.01` or marked pathogenic

##### 5.4 Target Region Variation
- Fetches all variants in phenotype-relevant genes

##### 5.5 VOI Phenotype Table
- Imports VOI phenotype data and filters by `VARID`

##### 5.6 Coverage
- Joins coverage data with target BED file
- Summarizes exon-level coverage and percent coverage

##### 5.7 Other Annotations
- **Sex prediction:** from `SEX` table
- **mtDNA haplogroup:** from `mtDNA` table
- **Ancestry:** from `ANCESTRY` table

##### 5.8 ClinVar Pathogenic Coverage
- Intersects ClinVar variants with coverage data
- Filters poorly covered or unconfirmed sites

##### 5.9 Founder Variants
- Joins with coverage and annotation if present


#### 6. Report Generation
- Uses `rmarkdown::render()` to create a sample-specific report:
  - Template: `bin/report.Rmd`
  - Output: `[SampleID_Flowcell]-[Phenotype/Description].report.html`
  - Output dir: `$WD/../REPORTS_<version>/<Flowcell>/`
  - Parameters: annotation tables, metadata, and config


#### 7. Final Output
- Writes processed samples to:  
  `Processed_SampleIDs.txt`




