# REBC_tools
Repo for code &amp; tasks used for REBC analysis

Terra tasks: 

* AllelicCapSeg_PP_CCF_fit_v3  calculates allelic copy number CCFs for each segment from an input AllelicCapSeg formated file. 

Inputs: 
    File ACS_FILE: output tsv file from task Model_Segments_PostProcessing_canonical_gatk414 
    String PURITY: tumor purity value 0<purity<=1
    String PLOIDY: tumor ploidy
    String ID: Tumor ID
Outputs:     
    File alleliccapseg_pp_ccf_fit_v3_out="${ID}.acs.ccf.tsv"
    File alleliccapseg_pp_ccf_fit_v3_collapsed_out="${ID}.acs.ccf.collapse.tsv"
