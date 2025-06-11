process ADDREADGROUP {
    tag "$sampleId - $flowcell"

      input:
    tuple val(sampleId), val(flowcell), path(bam_file) 

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.RG.bam"), emit: readgroup_bam_file

    script:
    """
    gatk AddOrReplaceReadGroups     \
        --INPUT ${bam_file}             \
        --OUTPUT "${sampleId}.RG.bam" \
        --RGID 4                    \
        --RGLB lib1                 \
        --RGPL illumina             \
        --RGPU unit1                \
        --RGSM ${sampleId}
    """

}
