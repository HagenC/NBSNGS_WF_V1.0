process BUILDSQL {
  tag "Building SQL database " 
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "SAMPLESforREPORTING.txt"
   
 
    input:
    val(sampleId)
   
    output:
    path('SAMPLESforREPORTING.txt'), emit: samples_reporting
    
    
 
    
 script:

    """
    buildsql.R $PWD  ${params.version}
    
     
    """
}
