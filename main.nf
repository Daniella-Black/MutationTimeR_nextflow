#! /usr/bin/env nextflow

//define channels from input file
//Channel 
//   .fromPath(params.inputlist)
//    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
//   .splitCsv(skip:1)
//    .map{sample, vcf_path,cnv_path, header, vcftobedpe -> [sample, vcf_path, cnv_path, header, vcftobedpe]}
 //   .set{ ch_input }


nextflow.enable.dsl=2  
//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample), file(vcf_path), file(cnv_path), file(header), file(vcftobedpe)

    output:
    file "*.tsv"
    file "*.txt"

    script:
    """
     CloudOS_MTR_input_script.R --vanilla $sample $vcf_path $cnv_path $header $vcftobedpe
    """
}

workflow {
    Channel.fromPath(params.inputlist) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample, file(row.vcf_path), file(row.cnv_path), file(row.header), file(row.vcftobedpe)) } 
        //| CloudOS_MTR_input
}