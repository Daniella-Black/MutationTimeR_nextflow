#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey,sv-> [tumour_sample_platekey, file(sv)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    //container = 'dockeraccountdani/mtrnf:latest' 
    container = 'dockeraccountdani/fitms2:latest' 
    //containerOptions '--volume ${workflow.workDir}/bin:/'
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'

    input:
    //set val(tumour_sample_platekey), file(somatic_cnv_vcf), file(somatic_small_variants_vcf_path) from ch_input
    set val(tumour_sample_platekey), file(sv)

    output:
    //file "small_variants_*.vcf.gz"
    //file "cnv_*.vcf.gz"
    file "sv_*.vcf.gz"
    //file "*_vaf_expected_vaf.pdf"
    //file "*_mT.csv"
    //file "*_mV.csv"
    //file "*_CLS.csv"
    //file "*.txt"

    script:
    """
    cp $sv sv_'$tumour_sample_platekey'.vcf.gz
    """ 
     //.map{tumour_sample_platekey,somatic_cnv_vcf,somatic_small_variants_vcf_path -> [tumour_sample_platekey, file(somatic_cnv_vcf), file(somatic_small_variants_vcf_path)]}
    //CloudOS_MTR_input_script.R '$tumour_sample_platekey' '$somatic_cnv_vcf' '$somatic_small_variants_vcf_path'
    //chmod +x $PWD/CloudOS_MTR_input_script.R
    //chmod +x bin/CloudOS_MTR_input_script.R
    //CloudOS_MTR_input_script.R
    //cp $somatic_small_variants_vcf_path small_variants_'$tumour_sample_platekey'.vcf.gz
    //cp $somatic_cnv_vcf cnv_'$tumour_sample_platekey'.vcf.gz
    // writeLines('$tumour_sample_platekey', paste0("out_", '$tumour_sample_platekey', ".txt"))
}
