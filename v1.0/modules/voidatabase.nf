process VOIDATABASE {
  tag "Generating VOIs database " 
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "VOI_DB_PHENOTYPE"
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "VOI_DB"
 
    input:
    val(sampleId)
   
    output:
    tuple path('VOI_DB'), path('VOI_DB_PHENOTYPE'), emit: VOI_databases
    
    
 
    
 script:

    """
    voidatabase.R $PWD 
    
     
    """
}
