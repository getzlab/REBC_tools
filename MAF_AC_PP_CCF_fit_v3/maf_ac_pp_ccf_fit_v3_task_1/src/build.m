pwd
rmpath	/home/unix/stewart/Projects/Cancer/tools/matlab
rmpath	/home/unix/stewart/Projects/tools/matlab
rmpath	/home/unix/stewart/Projects/matlab
rmpath	/home/unix/stewart/CancerGenomeAnalysis/trunk/matlab/mike
rmpath	/home/unix/stewart/CancerGenomeAnalysis/trunk/matlab/seq
rmpath	/home/unix/stewart/CancerGenomeAnalysis/trunk/matlab
addpath ../../../AllelicCapSeg_PP_CCF_fit_v3/alleliccapseg_pp_ccf_fit_v3_task_1/src
mcc -g -m fc_MAF_AC_PP_CCF_fit_v3;
quit;