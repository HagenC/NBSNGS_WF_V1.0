process MTDNA_HAPLOGROUP {
  tag "Predicting mtDNA haplogroup" 
   publishDir "${params.outputDir}/assets/temp/", mode: 'copy', overwrite: true, pattern: "mtDNA_haplotypedSamples.txt"
   
 
    input:
    val OPL_ch
   
    output:
    path('mtDNA_haplotypedSamples.txt'), emit: mtDNA_haplogroup
    
    
 
    
 script:

    """
    mtdnahaplogroup.R $PWD ${params.version}
    
     
    """
}
