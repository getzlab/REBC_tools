function [X1,X2]=fc_AllelicCapSeg_PP_CCF_fit_v3(AC,purity,ploidy,id,Gender)
%
%  AllelicCapSeg_PP_CCF_fit_v3 adds estimated tumor SCNA allele copy number [NA,NB] and
%  CCF's to AC based on puity and ploidy (from ABSOLUTE or other method).
%  This function makes estimates for all segments, even X and Y when gender
%  info is available from PAR.
%  AllelicCapSeg_PP_CCF_fit_v3 operates on one sample at a time (one purity, ploidy, gender)
%
%  AllelicCapSeg_PP_CCF_fit_v3 differs from v2 in that it fits CCFs to both alternate copy numbers.
%
%
%  Note: X output roughly correpsonds to the ABSOLUTE SEG_TAB file, with a few caveats
%      1) The model implements alternate (NA,NB) and normal allelic SCNA
%      state (1,1) for autosome where NA and NB can not be the normal copy
%      number.
%      2) NA and NB have corresponding CCF_NA and CCF_NB. A CCF of 0
%      implies the alternate allelic copy number does not exist in the
%      tumor sample
%      3) the normal CN occurs with fraction (1-CCF)
%      4) Can includes X and Y Chromosome estimates - making an assumption
%      of normal allelic copy number when there is no information from input AC
%      5) Can includes [NA NB] and [CCF_NA,CCF_NB] even for segments w/o het (mu_major/minor) info
%           when tau in NA segment within sigma_tau of neighbor(s)
%
%  inputs:
%     ALLELIC_CAPSEG_STRUCT fields:
%       tau, sigma_tau, mu_major, sigma_major, mu_minor,sigma_minor
%     purity: tumor sample fraction of cancer cells (0 to 1)
%     ploidy: cancer cell average CN over genome (0 to  inf)
%     gender: 'MALE' or 'FEMALE' ('M' or 'F' in any case combination)
%     PAR fields:
%       SIGMA_MIN: min sigma to be added in quadrature (0.01)
%       MAXCN: maximum tumor SCNA copy number (10)
%       SIGMA_THRESHOLD: signal:noise threshold for SCNA (2)
%
%  output:
%     X: AC struct with new fields:
%       cn_NA: SCNA Allele A copy number (lower or equal aCN)
%       cn_NB: SCNA Allele B copy number (higher or equal aCN)
%       cn_NA_CCF: CCF estimate NA
%       cn_NB_CCF: CCF estimate NA
%       cn_minor: fitted allele somatic copy number 
%       cn_major: fitted allele somatic copy number 
%       cn_minor_cluster: label for A allele somatic copy number cluster
%       cn_major_cluster: label for B allele somatic copy number cluster
%       IS_SCNA: tumor SCNA flag
%       purity:
%       ploidy:
%       gender: MALE or FEMALE
%       CNN: normal CN state (2, except 1 for male X or Y, 0 for female Y)
%
%  Chip 21 Mar 20
%

if ~isnumeric(purity)
   purity=str2double(purity)
end
if ~isnumeric(ploidy)
   ploidy=str2double(ploidy)
end

X1=AllelicCapSeg_PP_CCF_fit_v3(AC,purity,ploidy,Gender)
printStruct(X1,-1,[id '.acs.ccf.tsv'])

[X2]=collapse_AllelicCapSeg_PP_CCF_fit_v3(X1)
printStruct(X2,-1,[id '.acs.ccf.collapse.tsv'])

if isdeployed
    exit
end
