/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//Variant wrangling
include { VOIDATABASE                   } from '../modules/voidatabase.nf'
include { RENDER                        } from '../modules/render.nf'
include { BUILDSQL                      } from '../modules/buildsql.nf'



workflow REPORTING {
    take:
    wait_for_OPL_ch
    rmdFile




    main:
     //UPDATE_CLINVAR() 
     BUILDSQL(wait_for_OPL_ch.collect())

     BUILDSQL.out.samples_reporting
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row['SampleID_Flowcell'], row['Flowcell'], row['Description']) }
        .set { samples_for_reporting }

        //samples_for_reporting.view()
     VOIDATABASE(BUILDSQL.out.samples_reporting)
     RENDER(VOIDATABASE.out.VOI_databases)
     	
       

      //COMBINEVCF(samples_for_reporting, UPDATE_CLINVAR.out.clinvar_pathogenic_subset)  
      //VOIDATABASE(COMBINEVCF.out.wait.collect(), SAMPLESHEETS.out.samplesheet_collection)
      //PRERENDER(COMBINEVCF.out.combined_files, VOIDATABASE.out.voi_DB)
      //RENDER(PRERENDER.out.render_files, rmdFile)

}
