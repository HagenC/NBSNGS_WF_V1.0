/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



//variant calling
include { CLEANUP                       } from '../modules/cleanup.nf'
include { INDEX as INDEX_CLEAN          } from '../modules/index.nf'
include { INDELQUAL                     } from '../modules/indelqual.nf'
include { LOFREQNODEFAULT               } from '../modules/lofreq_nodefault.nf'
include { LOFREQCALL                    } from '../modules/lofreq_default.nf'
include { BQSR                          } from '../modules/bqsr.nf'
include { APPLY_BQSR                    } from '../modules/apply_bqsr.nf'
include { HAPLOTYPECALLER               } from '../modules/haplotypecaller.nf'
include { GATK_LA                       } from '../modules/leftAlign.nf'
include { INDEX as INDEX_HAPLOTYPECALLER} from '../modules/index.nf'
include { GTOVCF                        } from '../modules/gtovcf.nf'



//variant annotation
include { UPDATE_CLINVAR                } from '../modules/update_clinvar.nf'
include { ANNOTATE_CLINVAR              } from '../modules/annotate_clinvar.nf'
include { ANNOTATE_GNOMAD               } from '../modules/annotate_gnomad.nf'
include { ANNOTATE_ONTOLOGY             } from '../modules/annotate_ontology.nf'
include { OPL                           } from '../modules/opl.nf'

//QC
include { COVERAGE                      } from '../modules/coverage.nf'



workflow CALLING {
    take:
    bam_file_with_index




    main:


	CLEANUP(bam_file_with_index)
        INDELQUAL(CLEANUP.out.clean_bam_file)
        INDEX_CLEAN(INDELQUAL.out.bam_file)
        LOFREQNODEFAULT(INDEX_CLEAN.out.bam_file_w_index)
        LOFREQCALL(INDEX_CLEAN.out.bam_file_w_index)
        BQSR(CLEANUP.out.clean_bam_file)
        APPLY_BQSR(BQSR.out.bqsr_file)
        INDEX_HAPLOTYPECALLER(APPLY_BQSR.out.corrected_bam_file)
        COVERAGE(INDEX_HAPLOTYPECALLER.out.bam_file_w_index)
        HAPLOTYPECALLER(INDEX_HAPLOTYPECALLER.out.bam_file_w_index)
        //HAPLOTYPECALLER_VCF(INDEX_HAPLOTYPECALLER.out.bam_file_w_index)
        GTOVCF(HAPLOTYPECALLER.out.vcf_file)

   emit:
   wait_for_gvcf = GTOVCF.out.vcffile
}
