process MULTIQC {
    tag "MultiQC per flowcell"

    publishDir "${params.outputDir}/../QC_${params.version}/$flowcell/", pattern: "${flowcell}.multiqc_report.html", mode:'copy' 

    input:
     val(flowcell)
    path (fastq_metrics_files)
    path(markduplicate_metrics_file)

    output:
    path "${flowcell}.multiqc_report.html", emit: multiqc_html

    script:
    """
    multiqc $PWD/../QC_${params.version}/$flowcell --filename ${flowcell}.multiqc_report.html
    """

    
}
