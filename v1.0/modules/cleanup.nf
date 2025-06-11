process CLEANUP {
    tag "$sampleId - $flowcell"
    //publishDir "${params.outputDir}/SCheck/${flowcell}", mode: 'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file), path(bai)

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.clean.bam") , emit: clean_bam_file

    script:
    """
    samtools view                       \
        -F 1024                         \
        -F 512                          \
        -F 2048                         \
        -q 18                           \
        -b ${bam_file}                             \
        -o "${sampleId}.clean.bam"     \
       

    ## -F Filter based on the following flags:
    ## -F 1024: exclude read that are PCR or optical duplicate
    ## -F 512: exclude read that fails platform/vendor quality checks
    ## -F 2048: exclude supplementary alignments
    ## -q 20: exclude reads with mapping quality less than 20
    ## -b: output in BAM format
    """

 
}
