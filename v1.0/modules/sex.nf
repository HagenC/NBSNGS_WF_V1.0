process SEX {
  tag "Predicting Sex" 
  publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "SexProcessedIDs.txt"
  
  
  input:
    val OPL_ch
  
  output:
    path('SexProcessedIDs.txt'), emit: sex_pred
  
  
  
  
  script:
    
    """
    sexcheck.R $PWD  ${params.version}
    
     
    """
}
