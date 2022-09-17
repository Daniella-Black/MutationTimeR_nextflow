#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_small_variants_vcf_path,somatic_cnv_vcf -> [tumour_sample_platekey, file(somatic_small_variants_vcf_path), file(somatic_cnv_vcf)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    container = 'dockeraccountdani/mtrinputcloudos:latest' 
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_sample_platekey), file(somatic_small_variants_vcf_path), file(somatic_cnv_vcf) from ch_input

    output:
    //file "small_variants_*.vcf.gz"
    //file "cnv_*.vcf.gz"
    file "out_*.txt"

    script:
    """
    R  writeLines('$tumour_sample_platekey', paste0("out_", '$tumour_sample_platekey', ".txt"))
    """ 
    //cp $somatic_small_variants_vcf_path small_variants_'$tumour_sample_platekey'.vcf.gz
    //cp $somatic_cnv_vcf cnv_'$tumour_sample_platekey'.vcf.gz
}
