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
    tag"$tumour_sample_platekey"
    
    input:
    set val(tumour_sample_platekey), file(somatic_sv_vcf) from ch_input

    output:
    file "*.vcf.gz"
    
    script:
    """
    cp $somatic_sv_vcf somatic_sv_vcf_'$tumour_sample_platekey'.vcf.gz
    """ 
}
