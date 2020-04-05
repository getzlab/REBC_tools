function [X1,XA,A1,A2,AA]=fc_MAF_AC_PP_CCF_fit_v3(X,AC,purity,ploidy,id,Gender)
%
%  MAF_AC_PP_CCF_fit_v3 adds estimated mutation CCF and local tumor SCNA allele copy number
%  [NA,NB] and  CCF's to AC based on puity and ploidy (from ABSOLUTE or other method).
%  This function makes estimates for all segments, even X and Y when gender
%  info is available from PAR.
%  MAF_AC_PP_CCF_fit_v3 operates on one sample at a time (one purity, ploidy, gender)
%
%  Note: X output roughly correpsonds to the ABSOLUTE ABS_MAF file, with a
%  few caveats:
%      1) mutation multiplicity is 1 or NA or NB (no intermediate m) with a
%      prior that as the SCNA CCF goes up the likelihood that the SCNA preceeded
%      the mutation increases.
%      2) When allele fraction betapdf's exceed CCF=1, the overflow is assigned to CCF=1
%
%  inputs:
%     X: MAF struct with fields:
%       Chromosome, Start_position, t_alt_count, t_ref_count
%     T: ALLELIC_CAPSEG modified struct w/ fitted (NA,NB,NA_CCF,NB_CCF) fields:
%       tau, sigma_tau, mu_major, sigma_major, mu_minor,sigma_minor, NA, NB, NA_CCF, NB_CCF
%     purity: tumor sample fraction of cancer cells (0 to 1)
%     ploidy: cancer cell average CN over genome (0 to  inf)
%     gender: 'MALE' or 'FEMALE' ('M' or 'F' in any case combination)
%     PAR fields:
%       MAXCN: maximum tumor SCNA copy number (10)
%       SMOOTH_SPAN: taps in chisquare smoothing function (7)
%       DELTA: spacing in AF and CCF bins over range (0-1) and extended range for CCF
%       SIGMA_THRESHOLD: signal:noise threshold for SCNA (2)
%       CCF_PRIOR_SLOPE: factor favoring CCF=1 over CCF=0.5 linear
%
%  output:
%     X: MAF struct with new fields:
%       SCNA_NA: SCNA Allele A copy number (lower or equal CN)
%       SCNA_NB: SCNA Allele B copy number (higher CN)
%       SCNA_CCF_prob: vector of CCF probalities
%       SCNA_CCF_hat: CCF estimate (mode)
%       SCNA_CCF_Low: CCF 95% low CI
%       SCNA_CCF_High: CCF 95% high CI
%       IS_SCNA: tumor SCNA flag
%       purity:
%       ploidy:
%       gender: MALE or FEMALE
%       SCNA_NCN: normal CN state (2, except 1 for male X, Y, 0 for female Y)
%       dna_fraction_in_tumor: fraction of tumor cells in tumor sample
%       CCF_hat: CCF estimate
%       CCF_Low: CCF low CI (68%)
%       CCF_High: CCF high CI (68%)
%       CCF_CI95_low: CCF low CI (95%)
%       CCF_CI95_high: CCF high CI (95%)
%       CCF_mode: CCF mode
%       CCF_mean: average CCF
%       CCF_median: median CCF
%       i_tumor_f_corrected: alt allele fraction in tumor
%       pCCF: vector of CCF prob
%       clonal: clonal mutation flag
%
%
%  Chip 21 Mar 20
%

if ~isnumeric(purity)
   purity=str2double(purity)
end
if ~isnumeric(ploidy)
   ploidy=str2double(ploidy)
end

A1=AllelicCapSeg_PP_CCF_fit_v3(AC,purity,ploidy,Gender)
printStruct(A1,-1,[id '.acs.ccf.tsv'])

A2=collapse_AllelicCapSeg_PP_CCF_fit_v3(A1)
printStruct(A2,-1,[id '.acs.ccf.collapse.tsv'])

AA=AllelicCapSeg_PP_CCF_fit_v3_to_ABSOLUTE_SEGTAB(A2)
printStruct(AA,-1,'tmp.segtab.tsv')
cmd=[' sed ''1s/_Z_/\./g'' tmp.segtab.tsv > ' id '.AllelicCapSeg_PP_CCF_fit_v3.ABSOLUTE.segtab.txt']
unix(cmd)


[ X1, XA] = MAF_AC_PP_CCF_fit_v3(X,A2,purity,ploidy,Gender)
printStruct(X1,-1,[id '.MAF_AC_PP_CCF_fit_v3.maf'])
printStruct(XA,-1,'tmp.ABSOLUTE.maf')
cmd=[' sed  ''1s/_Z_/\./g'' tmp.ABSOLUTE.maf | sed ''1s/pccf_0Z/0\./g'' | sed  ''1s/pccf_0/0/g'' | sed ''1s/pccf_1/1/g'' > ' id '.MAF_AC_PP_CCF_fit_v3.ABSOLUTE.maf ']
unix(cmd)
% cmd=[' sed -i ''1s/pccf_0Z/0\./g'' tmp.ABSOLUTE.maf ']
% unix(cmd)
% cmd=[' sed -i ''1s/pccf_0/0/g'' tmp.ABSOLUTE.maf ']
% unix(cmd)
% cmd=[' sed ''1s/pccf_1/1/g'' tmp.ABSOLUTE.maf > ' id '.MAF_AC_PP_CCF_fit_v3.ABSOLUTE.maf']
% unix(cmd)


if isdeployed
    exit
end

end
%%
function test
clear
rmpath('/Users/stewart/CancerGenomeAnalysis/trunk/matlab')
rmpath('/Users/stewart/CancerGenomeAnalysis/trunk/matlab/seq')
rmpath('/Volumes/GoogleDrive/My Drive/Cancer/tools/matlab')
rmpath('/Users/stewart/Projects/Cancer/tools/matlab')
rmpath('/Users/stewart/Projects/matlab')
addpath('/Users/stewart/Projects/REBC_tools/MAF_AC_PP_CCF_fit_v3/maf_ac_pp_ccf_fit_v3_task_1/src/')
addpath('/Users/stewart/Projects/REBC_tools/AllelicCapSeg_PP_CCF_fit_v3/alleliccapseg_pp_ccf_fit_v3_task_1/src/')

    
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/Mutations/test')
pid1='REBC-ACAB-TP-NB'
%pid1='REBC-ACAE-TP-NT'
X1=load_tsv([pid1 '_consensus_fix_NALT01.maf'])
A=load_tsv('../../CN/REBC_primary_pair_393.cnv_postprocessing_tumor_acs.tsv')
P=load_tsv('../../FC/REBC.pairs.15Mar2020.tsv')
PT=load_tsv('../../FC/REBC.participants.23Feb2020.tsv')
A1=trimStruct(A,ismember(A.sample,pid1))
P1=trimStruct(P,ismember(P.entity_pair_id,pid1))
PT1=trimStruct(PT,ismember(PT.entity_participant_id,P1.participant))
purity=P1.purity(1)
ploidy=P1.absolute_ploidy(1)
gender=PT1.Gender;
[X1,XA,A1,A2,AA] = fc_MAF_AC_PP_CCF_fit_v3(X1,A1,purity,ploidy,pid1,gender)


% reset path
startup
end

