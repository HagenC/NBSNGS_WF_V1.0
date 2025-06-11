process FASTQC {
    tag "$flowcell"
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/fastqc/", mode: 'copy'

    input:
    input:
    tuple val(sampleId), path(fastq1), path(fastq2), val(flowcell)
    
    

    output:
    path "*.zip", emit: zipFiles
    path "*.html", emit: htmlReports

    script:
    """
    fastqc --threads ${task.cpus}  ${fastq1} ${fastq2}
    """
}


