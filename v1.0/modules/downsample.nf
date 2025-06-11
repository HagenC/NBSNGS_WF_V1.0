process DOWNSAMPLE {
    tag "$sample_Id"
     //publishDir "${params.outputDir}/SCheck/${project}", pattern: "${sample_Id}_R1.fq.gz", mode:'copy'
     //publishDir "${params.outputDir}/SCheck/${project}", pattern: "${sample_Id}_R2.fq.gz", mode:'copy'
    
    
    input:
    tuple val(sample_Id), path(fastq_1) , path(fastq_2), val(project)
    

    output:
    tuple val(sample_Id), path("${sample_Id}_R1.fq.gz"), path("${sample_Id}_R2.fq.gz"), val(project), emit: sampled_reads

    script:
    
    """
    # 100× coverage × 34,156,490 bp ÷ 100 bp/read = 34,156,490 reads
    reformat.sh \
        in1=${fastq_1} \
        in2=${fastq_2} \
        out1="${sample_Id}_R1.fq.gz" \
        out2="${sample_Id}_R2.fq.gz" \
        sampleseed=1 \
        samplereadstarget=93487190\
    """
   

}

