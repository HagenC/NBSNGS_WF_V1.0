#!/usr/bin/env nextflow

process ALIGNMENT {
    tag "${sampleId} - $flowcell"
   //publishDir "${params.outputDir}/SCheck/${flowcell}",  pattern: "${sampleId}.sam", mode:'copy'
   publishDir "${params.outputDir}/../log_${params.version}/${flowcell}/mapping/", pattern: "${sampleId}.log", mode:'copy'

    input:
    tuple val(sampleId), path(fastq1), path(fastq2), val(flowcell)

    output:
    tuple val(sampleId),val(flowcell), path("${sampleId}.bam"), emit: raw_bam_file
    path "${sampleId}.log"


    script:
    """
     bwa-mem2 mem  \
        -t ${task.cpus}                                \
        -K 1000000000                                  \
       ${params.reference_genome}   \
        $fastq1 $fastq2 2> >(tee -a "${sampleId}.log" >&2) | samtools sort -@ ${task.cpus} -o - - | samtools view -L ${params.call_bed} -b -o ${sampleId}.bam -
	

    """
}



//2> >(tee -a "${sampleId}.log" >&2) | samtools sort -@ ${task.cpus} -o - - | samtools view -L ${params.bedFile_padded} -b -o ${sampleId}.bam -
