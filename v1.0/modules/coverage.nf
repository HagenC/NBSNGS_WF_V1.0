process COVERAGE {
    tag "Coverage - $sampleId  - $flowcell"
    publishDir "${params.outputDir}/../QC_${params.version}/${flowcell}/", pattern: "${sampleId}.cov.gz", mode: 'copy'
    
    
    input:
    tuple val(sampleId), val(flowcell), path(bam_file), path(bam_bai)

    output:
    path("${sampleId}.cov.gz")

    script:
    """
    samtools depth -a -q 18 -Q 20 -b ${params.coverage_bed} $bam_file > "${sampleId}.cov" 
    gzip -f "${sampleId}.cov"
       
    """
}



