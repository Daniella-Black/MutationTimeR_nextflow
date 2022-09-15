#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{sample, vcf_path,cnv_path, header, vcftobedpe -> [sample, file(vcf_path), file(cnv_path), file(header), file(vcftobedpe)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set sample, file(vcf_path), file(cnv_path), file(header), file(vcftobedpe) from ch_input

    output:
    //file "*.tsv"
    //file "*.txt"
    file "x_*"

    script:
    """
    x_$sample <- gunzip $vcf_path
    """ 
}
