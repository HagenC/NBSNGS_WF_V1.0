process ANCESTRY {
  tag "Predicting Ancestry" 
  publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "ancestry_samples.txt"
  
  
  input:
    val OPL_ch
  
  output:
    path('ancestry_samples.txt'), emit: ancestry_pred
  
  
  
  
  script:
    
    """
    ancestrymodel.R $PWD ${params.ancestry_reference} ${params.ancestry_RS_selection} ${params.version}
    
     
    """
}
