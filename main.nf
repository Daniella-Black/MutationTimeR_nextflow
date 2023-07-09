#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey,somatic_small_variants_vcf_path, somatic_cnv_vcf,tumour_purity, header-> [tumour_sample_platekey, file(somatic_small_variants_vcf_path), file(somatic_cnv_vcf), tumour_purity, file(header)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    maxForks 900
    //errorStrategy 'ignore'
    maxRetries 3
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_sample_platekey), file(somatic_small_variants_vcf_path), file(somatic_cnv_vcf), val(tumour_purity), file(header) from ch_input

    output:
    file "*_mutationtimer_input_CNVs.txt"
    file "*_mutationtimer_input_SNVs.txt"


    script:
    """
    CloudOS_MTR_input_script.R '$tumour_sample_platekey' '$somatic_small_variants_vcf_path' '$somatic_cnv_vcf' '$tumour_purity' '$header'
    """ 
}
