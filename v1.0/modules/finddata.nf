process FINDDATA {
  tag "Collecting and matching fastq-files" 
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "Sample.collection.info"
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "ProcessedSamples.txt" 
   //publishDir "${params.outputDir}/assets/sanity/", mode: 'copy', overwrite: true, pattern: "SnpEff_geneFilter.list_observation_count.tsv"  
   
  
  output:
    path 'Sample.collection.info'
    path 'ProcessedSamples.txt'
    path 'ProccessingIDs.txt', emit: alignment_tuple
    //path 'SnpEff_geneFilter.list_observation_count.tsv'
   
    
 script:

    """
    module load R/4.3.1
    sqlsamplesheet.R $PWD ${params.coverage_bed} ${params.call_bed} ${params.target_bed} ${params.snpeff_gtf} ${params.version}
     
    """
}


