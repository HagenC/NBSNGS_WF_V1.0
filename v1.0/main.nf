#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """\
    ========================================================================================
             NBSNGS WF  v.1.0  
    =======================================================================================
    
    Regenerate alignment, QC? : REMEMBER to delete blocker-files (*.GATK.g.vcf.gz)
    
    Regenerate reports? : REMEMBER to delete blocker-files (*.html)
    
    See STRUCTURE.ME for structure.
    
    Debugging?: set conditional to false.
    
    """
    .stripIndent()

/*

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONDITIONAL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


params.delete_work_dir = true  // true or false

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//mapping 
include { COLLECTDATA          	} from './subworkflows/collectdata.nf'
include { MAPPING          	} from './subworkflows/mapping.nf'
include { CALLING          	} from './subworkflows/calling.nf'
include { ANNOTATING          	} from './subworkflows/annotating.nf'
include { REPORTING          	} from './subworkflows/reporting.nf'
//include { UPDATE_CLINVAR        } from './modules/update_clinvar.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LINKS and PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


params.outputDir = "."
rmdFile = "${PWD}/bin/report.Rmd"




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//main workflow
workflow  {
    NbsNgsWF()
   
    
    
   
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NbsNgsWF {
        SENDMAIL()
	//UPDATE_CLINVAR()
	ALIGNMENT = COLLECTDATA().alignment_tuple
        
	
	
	ALIGNMENT
    	.splitCsv(header: true, sep: '\t')
    	.map { row -> tuple(row['SampleID_Flowcell'], file(row['fastq_1']), file(row['fastq_2']), row['Flowcell']) }
    	.set { samples_for_alignment }
        //samples_for_alignment.view()
	
	samples_for_alignment
        .map { tuple -> tuple[3] }  
        .distinct()                 
        .set { unique_flowcells }  
	 //unique_flowcells.view() 
	
	 	
	
	bam_file_w_index_ch = MAPPING(samples_for_alignment, unique_flowcells).bam_file_with_index
	
	

	if(params.doCalling) {
    	wait_for_gvcf_ch = CALLING(bam_file_w_index_ch).wait_for_gvcf
 }

	
	 if(params.doReporting){
   wait_for_ANCESTRY_ch = ANNOTATING(wait_for_gvcf_ch).wait_for_ANCESTRY
   .ifEmpty { Channel.of("Continue anyway") }
  REPORTING(wait_for_ANCESTRY_ch, rmdFile)
       }
	 

	
      

}


workflow.onComplete {
    
  
    
    if (workflow.success) {
        log.info "\nNbsNgs WF has successfully completed!\n"
	
       
	
	
//MAIL:
        def msg = """\
        NbsNgs WF execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
         
        sendMail(subject: "NBSNGS WF v.1.0 sucecesfully completed", body:  msg ,  attach: "./assets/temp/ProcessedSamples.txt", to: 'unc@ssi.dk')

//Delete workdir:

        if (params.delete_work_dir) {
            def deleteCmd = "rm -rf ${params.outputDir}/work"
            def process = deleteCmd.execute()
            process.waitFor()
            if (process.exitValue() == 0) {
                log.info "Work directory deleted based on true-setting!."
            } else {
                log.error "Failed to delete work directory."
            }
        } else {
            log.info "Work directory NOT executed based on true-setting!."
        }
//ON failure:
    } else {
        log.error "Damn! .. something went wrong"
        sendMail(subject: "NBSNGS WF v.1.0 - FAILED", body: "NBSNGS WF v.1.0 - FAILED!", attach: "./assets/temp/ProcessedSamples.txt", to: 'unc@ssi.dk')
    }
}


process SENDMAIL {

  

    output:
    path "email_notification.txt"

    script:
    """
    echo "NBSNGS WF v.1.0 started" | mail -s "NBSNGS WF v.1.0  started" , unc@ssi.dk
    echo "Email sent to unc@ssi.dk" > email_notification.txt
    """
}





