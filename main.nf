#! /usr/bin/env nextflow

//define channels from input file
//Channel 
//   .fromPath(params.inputlist)
//    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
//   .splitCsv(skip:1)
//    .map{sample, vcf_path,cnv_path, header, vcftobedpe -> [sample, vcf_path, cnv_path, header, vcftobedpe]}
 //   .set{ ch_input }



workflow {
    Channel.fromPath(params.inputlist) \
        | splitCsv(header:true) \
        | map { row-> tuple(row.sample, file(row.vcf_path), file(row.cnv_path), file(row.header), file(row.vcftobedpe)) } 
        //| CloudOS_MTR_input
}
