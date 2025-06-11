process INDELQUAL {
    tag "$sampleId - $flowcell"
     //conda "/home/dksund.dk/unc/miniconda3/envs/lofreq_env"
    
    //publishDir "${params.outputDir}/../log_${params.version}/lofreq_indelqual/${flowcell}/", pattern: "${sampleId}.lofreq.iq.log", mode:'copy'

    input:
    tuple val(sampleId), val(flowcell), path(bam_file)
   

    output:
    tuple val(sampleId), val(flowcell), path("${sampleId}.iq.bam"),  emit: bam_file
    path "${sampleId}.lofreq.iq.log"

    script:
    """
    lofreq indelqual                     \
        --dindel                                        \
        -f ${params.reference_genome}                   \
        -o "${sampleId}.iq.bam"                         \
	 ${bam_file}                                      \
	 > >(tee -a "${sampleId}.lofreq.iq.log")
    """
}
