#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey,somatic_sv_vcf -> [tumour_sample_platekey, file(somatic_sv_vcf)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    container = 'dockeraccountdani/fitms2:latest' 
    //containerOptions '--volume ${workflow.workDir}/bin:/'
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}/$tumour_sample_platekey", mode: 'copy'

    input:
    set val(tumour_sample_platekey), file(somatic_sv_vcf) from ch_input

    output:
    //file "small_variants_*.vcf.gz"
    file "sv_*.vcf.gz"
    //file "*_vaf_expected_vaf.pdf"
    //file "*_mT.csv"
    //file "*_mV.csv"
    //file "*_CLS.csv"
    //file "*.txt"

    script:
    """
    cp $somatic_sv_vcf '$tumour_sample_platekey'.vcf.gz
    """ 
    //chmod +x $PWD/CloudOS_MTR_input_script.R
    //chmod +x bin/CloudOS_MTR_input_script.R
    //CloudOS_MTR_input_script.R
    //cp $somatic_small_variants_vcf_path small_variants_'$tumour_sample_platekey'.vcf.gz
    //cp $somatic_cnv_vcf cnv_'$tumour_sample_platekey'.vcf.gz
    // writeLines('$tumour_sample_platekey', paste0("out_", '$tumour_sample_platekey', ".txt"))
}
