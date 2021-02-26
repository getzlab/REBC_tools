# REBC_tools
Repo for code &amp; tasks used for REBC analysis

Terra tasks: Each with matlab src directory, Dockerfile, and workflow definition  WDL for Terra. 

* AllelicCapSeg_PP_CCF_fit_v3  calculates allelic copy number CCFs for each segment from an input AllelicCapSeg formated file. Additionally the tool provides a second output is a "smoothed" version of the allelic copy number CCF file that was  

Inputs: 
    File ACS_FILE: output tsv file from task Model_Segments_PostProcessing_canonical_gatk414 
    String PURITY: tumor purity value 0<purity<=1
    String PLOIDY: tumor ploidy
    String ID: Tumor ID
    
Outputs:     
    File alleliccapseg_pp_ccf_fit_v3_out:  includes fields with allelic copy number CCF estimates 
    File alleliccapseg_pp_ccf_fit_v3_collapsed_out: same output format as above,  
