#! /usr/bin/env nextflow
//define channels from input file
Channel 
    .fromPath(params.inputlist)
    .ifEmpty {exit 1, "Cannot find input file : ${params.inputlist}"}
    .splitCsv(skip:1)
    .map{tumour_sample_platekey, somatic_small_variants_vcf_path.x,somatic_cnv_vcf -> [sample, file(somatic_small_variants_vcf_path.x), file(somatic_cnv_vcf)]}
    .set{ ch_input }


//run the script to make MTR input on above file paths
process  CloudOS_MTR_input{
    tag"$sample"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    set sample, file(somatic_small_variants_vcf_path.x), file(somatic_cnv_vcf) from ch_input

    output:
    file "small_variants_*.vcf.gz"
    file "cnv_*.vcf.gz"

    script:
    """
    cp $somatic_small_variants_vcf_path.x small_variants_'$sample'.vcf.gz
    cp $somatic_cnv_vcf cnv_'$sample'.vcf.gz
    """ 
}
