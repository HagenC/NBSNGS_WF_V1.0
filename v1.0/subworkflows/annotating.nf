/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//variant annotation
include { GATK_LA                       } from '../modules/leftAlign.nf'
//include { UPDATE_CLINVAR                } from '../modules/update_clinvar.nf'
include { ANNOTATE_CLINVAR              } from '../modules/annotate_clinvar.nf'
include { ANNOTATE_GNOMAD               } from '../modules/annotate_gnomad.nf'
include { ANNOTATE_ONTOLOGY             } from '../modules/annotate_ontology.nf'
include { OPL                           } from '../modules/opl.nf'
include { UPDATE_CLINVAR                } from '../modules/update_clinvar.nf'
include { MTDNA_HAPLOGROUP              } from '../modules/mtdna_haplogroup.nf'
include { SEX                           } from '../modules/sex.nf'
include { ANCESTRY                      } from '../modules/ancestrymodel.nf'

workflow ANNOTATING {
    
    
    
    take:
    gtovcf_file


   

    main:
      UPDATE_CLINVAR() 
      GATK_LA(gtovcf_file) 
      ANNOTATE_GNOMAD(GATK_LA.out.la_vcf_file) 
      ANNOTATE_ONTOLOGY(ANNOTATE_GNOMAD.out.CLNGNOMAD_la_vcf_file) 
      ANNOTATE_CLINVAR(ANNOTATE_ONTOLOGY.out.CLNGNOMADSO_la_vcf_file, UPDATE_CLINVAR.out.clinvar_vcf, UPDATE_CLINVAR.out.clinvar_tbi) 
      OPL(ANNOTATE_CLINVAR.out.CLN_la_vcf_file) 
      OPL_ch = OPL.out.OPL_vcf_file.collect() 
      MTDNA_HAPLOGROUP(OPL_ch) 
      SEX(OPL_ch) 
      ANCESTRY(OPL_ch) 

   emit:
   wait_for_ANCESTRY = ANCESTRY.out.ancestry_pred
}
