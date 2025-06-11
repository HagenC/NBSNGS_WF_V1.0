process RENDER {
  tag "Render reports" 
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "Processed_SampleIDs.txt"

 
   input:
   tuple path(VOI_DB), path(VOI_DB_PHENOTYPE)
   
   output:
   path("Processed_SampleIDs.txt")
   
    
 
    
 script:

    """
    mkdir -p ../REPORTS_${params.version}
    render.R $PWD $VOI_DB $VOI_DB_PHENOTYPE ${params.version}
     
    """
}
