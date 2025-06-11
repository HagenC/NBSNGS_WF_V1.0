/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAPPING SUB-WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                    } from '../modules/fastqc'
include { ALIGNMENT                 } from '../modules/alignment'
include { INDEX  as RAW_INDEX       } from '../modules/index'
include { MOSDEPTH  as RAW_DEPTH    } from '../modules/mosdepth'
include { ADDREADGROUP              } from '../modules/addreadgroup'
include { HSMETRICS                 } from '../modules/hs_metrics'
include { MARKDUPLICATES            } from '../modules/markduplicates'
include { ALIGNMENT_METRICS         } from '../modules/alignment_metrics'
include { INSERT_SIZE_METRICS       } from '../modules/insert_size_metrics'
include { MULTIQC                   } from '../modules/multiqc'
include { INDEX                     } from '../modules/index'
include { DOWNSAMPLE                } from '../modules/downsample'

workflow MAPPING {
    take:
    sampletable
    flowcells




    main:

        DOWNSAMPLE(sampletable)
        FASTQC(DOWNSAMPLE.out.sampled_reads)
        qc_ch = FASTQC.out.zipFiles

     
    ALIGNMENT(DOWNSAMPLE.out.sampled_reads)
    RAW_INDEX(ALIGNMENT.out.raw_bam_file)
    ADDREADGROUP(ALIGNMENT.out.raw_bam_file)
    MARKDUPLICATES(ADDREADGROUP.out.readgroup_bam_file)
    INDEX(MARKDUPLICATES.out.marked_bam_file)

    if(params.doQC){
    RAW_DEPTH(RAW_INDEX.out.bam_file_w_index)
    HSMETRICS(ALIGNMENT.out.raw_bam_file)
    ALIGNMENT_METRICS(ALIGNMENT.out.raw_bam_file)
    //GC_METRICS(ALIGNMENT.out.raw_bam_file, "")
    INSERT_SIZE_METRICS(ALIGNMENT.out.raw_bam_file)
    QC_ch = qc_ch.mix(FASTQC.out.htmlReports, 
                     HSMETRICS.out.metrics_file, ALIGNMENT_METRICS.out.metrics_file,INSERT_SIZE_METRICS.out.metrics_file, 
		     RAW_DEPTH.out.per_base_depth, RAW_DEPTH.out.region_depth, RAW_DEPTH.out.global_dist, RAW_DEPTH.out.region_dist, RAW_DEPTH.out.summary ).collect()
   MULTIQC(flowcells, QC_ch, MARKDUPLICATES.out.metrics_file)
    }else{
   MULTIQC(flowcells, qc_ch, MARKDUPLICATES.out.metrics_file)

  }
   emit:
   bam_file_with_index = INDEX.out.bam_file_w_index
}
