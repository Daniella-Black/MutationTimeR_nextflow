#! /usr/bin/env nextflow

//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{sample, vcf_path,cnv_path, header, vcftobedpe -> [sample, vcf_path, cnv_path, header, vcftobedpe]}
    .set{ ch_input }

//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set val(sample), val(vcf_path), val(cnv_path), val(header), val(vcftobedpe) from ch_input

    output:
    file "my_data.csv"
    //file "*.tsv"
    //file "*.txt"
    

    script:
    """
    CloudOS_MTR_input_script.R $sample $vcf_path $cnv_path $header $vcftobedpe 
    """
    
    //output:
    //path '*.txt'
    
    //script:
    //'''
    //echo \$sample > file1.txt
    //echo \$vcf_path > file2.txt
    //echo \$cnv_path > file3.txt
    //'''
}
