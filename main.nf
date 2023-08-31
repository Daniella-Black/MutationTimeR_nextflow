#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_cnv_vcf,tumour_purity,organ-> [tumour_sample_platekey, file(somatic_cnv_vcf), tumour_purity, organ]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    maxForks 900
    errorStrategy 'ignore'
    maxRetries 3
    tag"$tumour_sample_platekey"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(tumour_sample_platekey), file(somatic_cnv_vcf), val(tumour_purity), val(organ) from ch_input

    output:
    file "*_CNVs.tsv"
    //file "*_SNVs.txt"
    //file '*_unique_filter_fields.tsv'


    script:
    """
    CloudOS_MTR_input_script.R '$tumour_sample_platekey' '$somatic_cnv_vcf' '$tumour_purity' '$organ'
    """ 
}
