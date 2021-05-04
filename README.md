# REBC_tools
Repo for code &amp; tasks used for REBC analysis


Terra tasks: Each with matlab src directory, Dockerfile, and workflow definition  WDL for Terra. 
---

* AllelicCapSeg_PP_CCF_fit_v3  calculates allelic copy number CCFs for each segment from an input AllelicCapSeg formated file. Additionally the tool provides a second output is a "smoothed" version of the allelic copy number CCF file that was  

  * Inputs: 
    * File ACS_FILE: output tsv file from task Model_Segments_PostProcessing_canonical_gatk414 
    * String PURITY: tumor purity value 0<purity<=1
    * String PLOIDY: tumor ploidy
    * String ID: Tumor ID
    
  * Outputs:     
    * File alleliccapseg_pp_ccf_fit_v3_out:  includes fields with allelic copy number CCF estimates 
    * File alleliccapseg_pp_ccf_fit_v3_collapsed_out: smoothed segments in same output format as above.  

* MAF_AC_PP_CCF_fit_v3  calculates mutation CCFs from an input MAF and allelicCapSeg formated file. 

  * Inputs: 
    * File MAF_FILE: output MAF file from task filter_NALT01_maf
    * File ACS_FILE: output tsv file from task Model_Segments_PostProcessing_canonical_gatk414 
    * String PURITY: tumor purity value 0<purity<=1
    * String PLOIDY: tumor ploidy
    * String Gender: male or female for X chromosome CCF estimate
    * String ID: Tumor ID
    
  * Outputs:     
    * File alleliccapseg_pp_ccf_fit_v3_out:  includes fields with allelic copy number CCF estimates 
    * File alleliccapseg_pp_ccf_fit_v3_collapsed_out: smoothed segments in same output format as above.  
    * File alleliccapseg_pp_ccf_fit_v3_ABSOLUTE_SEGTAB: allelic copy number file in ABSOLUTE format 
    * File maf_ac_pp_ccf_fit_v3_maf:  MAF file with fields related to mutation CCF estimates and distributions. 
    * File maf_ac_pp_ccf_fit_v3_ABSOLUTE_maf:  MAF file with fields related to mutation CCF estimates in ABSOLUTE format. 

* SV_CCF_v3  calculates SV CCFs from an input SV tsv in breakpointer format and allelicCapSeg formated file. 

  * Inputs: 
    * File SV_BP_FILE: output BP tsv file from task REBC_SV_consensus_filter
    * File ACS_FILE: output tsv file from task Model_Segments_PostProcessing_canonical_gatk414 
    * String PURITY: tumor purity value 0<purity<=1
    * String PLOIDY: tumor ploidy
    * String Gender: male or female for X chromosome CCF estimate
    * String ID: Tumor ID
    
  * Outputs:     
    * File SV_BP_CCF_v3_tsv:  BP tsv includes fields with SV CCF estimates 


Matlab scripts for post-processing outside of Terra: Most of these are matlab utility functions used by a main analysis script. 
---

* alleleFraction_clonal_hets_Purity_estimate.m: makes purity estimate based on finding somatic clonal het peak. 
* REBC_SV_hotspots_4May2020.m:  Mark or remove SV events with exactly recurrent breakpoints - likely to be mapping artifacts. 
* REBC_balanced_inv2_clusters_unfiltered_27Jun2020.m: clusters SV breakpoints into balanced, inversion, and complex chains of SVs. 
* revise_GISTIC_sample_tables.m: Correct GISTIC2 output tables marking samples for arm-level and focal SCNAs. 
* REBC_supplemental_figures_mutations.m: matlab script that generated supplmental figures 7 and 14. 
    
Optional python script to identify dupicated SV events - not actually needed for SV processing.
---
*  dedup_svs_v3.sh, which calls dedup_svs_v3.py
