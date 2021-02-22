
% gsutil -m cp gs://fc-035f5652-acf7-4642-abb7-e8c10848c8ed/f2ddbf3d-6a23-4713-8f16-51f121dc41dc/maf_aggregator_workflow/de4b80d0-6f1c-4383-a032-edebf4fe1925/call-maf_aggregator_task/REBC_primary_pair_393.aggregated.maf  /tmp/.

%  cut -f1,5-6,9-10,15-16,42,45,46,53,89,90,102,108,114,135,137,139,141,143,150,171,180,183,193,196,199,201  /tmp/REBC_primary_pair_393.aggregated.maf > /tmp/REBC_primary_pair_393.all.trim.maf

clear
cd ~/GoogleDrive/Cancer/REBC
X=load_tsv('Mutations/REBC_primary_pair_393.all.trim.maf')

T=load_tsv('MainTable/REBC_THCA_TABLE_V16.27Jul2020.tsv')
T=trimStruct(T,T.WGS_TP_N_QCPASS>0)
T=trimStruct(T,strfindk(T.LABELED_REBC_ID,'REBC')>0)
X.LABELED_REBC_ID=X.Tumor_Sample_Barcode;
X.LABELED_REBC_ID=regexprep(X.LABELED_REBC_ID,'-YQ-','-');
X.LABELED_REBC_ID=cellfun(@(x) x(1),regexp(X.LABELED_REBC_ID,'-TTP','split'));
X.LABELED_REBC_ID=cellfun(@(x) x(1),regexp(X.LABELED_REBC_ID,'-TP','split'));

% only analyze pass QC set 383
X=trimStruct(X,ismember(X.LABELED_REBC_ID,T.LABELED_REBC_ID))

tab(X.LABELED_REBC_ID)
% 
% myhist(X.t_alt_count,0:50)
% myhist(X.t_alt_count(ismember(X.Variant_Type,{'INS','DEL'})),0:50)  
% tab(X.t_alt_count<1)  % 1248 Strelka 0's 
% printStruct(X,X.t_alt_count<1)


X=trimStruct(X,X.t_ref_count>=0)


X.vaf=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
[i m]=ismember(X.LABELED_REBC_ID,T.LABELED_REBC_ID);
X.contEst1(i,1)=T.WGS_TP_N_contEst_contamination(m(i))/100.0;
X.pcontEst1=betacdf(X.contEst1,X.t_alt_count+1,X.t_ref_count+1)

% k1=find(abs(X.contEst1-0.001)<1e-5);
% k2=find(abs(X.contEst1-0.002)<1e-5);
% k3=find(abs(X.contEst1-0.003)<1e-5);
% plot(X.vaf(k1),X.pcontEst1(k1),'.',X.vaf(k2),X.pcontEst1(k2),'.',X.vaf(k3),X.pcontEst1(k3),'.')

PTHRESHOLD=1e-5
X.contEst_PASS=X.pcontEst1<PTHRESHOLD;

tab(X.contEst_PASS)

tab(X.realign_judgment)
X.blat_PASS=ismember(X.realign_judgment,'KEEP');
X.oxoG_PASS=1*(X.i_oxoG_cut<1);
X.M1=1*(X.M1>0)
X.M2=1*(X.M2>0)
X.SvABA=1*(contains(X.SvABA,'1'))
X.Snowman=1*(contains(X.Snowman,'1'))
X.Strelka1=1*(contains(X.Strelka1,'1'))
X.Strelka2=1*(contains(X.Strelka2,'1'))
X.vote=X.M1+X.M2+1*(X.Strelka1|X.Strelka2)+1*(X.SvABA|X.Snowman)

X.consensus_PASS=1*(X.vote>1);
X.NALT_PASS=1*(X.n_alt_count<=1);

tab(X.NALG)

myhist(X.n_alt_count,0:10)


figure(1)
clf
ALGS={'M1','M2','Strelka1','Strelka2','Snowman','SvABA'}
P=[]
P.log=1
P.angle=25
[h,sv] =  upset_plot(X,ALGS,P);
set(h(1),'ylim',[0.9 5e6])
AREA='Doc/draft-20201009/plots/'
PSET='REBC_primary_pair_383'
f=[AREA PSET '.Mut_raw_algs_UpSet.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')

figure(2)
passfilt=(X.pcawg_pon_pass_loglike>0)&(X.rebc_pon_pass_loglike>0)&(X.blat_PASS>0)&(X.oxoG_PASS>0)&(X.NALT_PASS>0);
tab(passfilt)
XF=trimStruct(X,passfilt)
clf
ALGS={'M1','M2','Strelka1','Strelka2','Snowman','SvABA'}
P=[]
P.log=1
P.angle=25
[h,sv] =  upset_plot(XF,ALGS,P);
set(h(1),'ylim',[0.9 5e6])
AREA='Doc/draft-20201009/plots/'
PSET='REBC_primary_pair_383'
f=[AREA PSET '.Mut_filt_algs_UpSet.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')


X.PON_REJECT=(X.pcawg_pon_pass_loglike.*X.rebc_pon_pass_loglike)<1;
X.BLAT_REJECT=X.blat_PASS<1;
X.CONTEST_REJECT=X.contEst_PASS<1;
X.CONSENSUS_REJECT=X.consensus_PASS<1;
X.NALT_REJECT=X.NALT_PASS<1;
X.OXOG_REJECT=zeros(size(X.oxoG_PASS));

X.PASS=1*((X.pcawg_pon_pass_loglike.*X.rebc_pon_pass_loglike.*X.blat_PASS.*X.contEst_PASS.*X.consensus_PASS.*X.NALT_PASS)>0);
tab(X.PASS)


tab(isnan(X.pcawg_pon_pass_loglike))
tab(isnan(X.rebc_pon_pass_loglike))
tab(isnan(X.blat_PASS))
tab(isnan(X.contEst_PASS))
tab(isnan(X.consensus_PASS))
tab(isnan(X.NALT_PASS))
tab(isnan(X.PASS))


FILT={'PASS','PON_REJECT','BLAT_REJECT','CONTEST_REJECT','CONSENSUS_REJECT','NALT_REJECT','OXOG_REJECT'}
figure(3)
clf
P=[]
P.log=1
P.angle=25
[h,sv] =  upset_plot(X,FILT,P);
set(h(1),'ylim',[0.9 5e6])
f=[AREA PSET '.Mut_raw_filters_UpSet.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')


figure(4)
passfilt=(X.consensus_PASS>0)&(X.pcawg_pon_pass_loglike>0)&(X.rebc_pon_pass_loglike>0)&(X.blat_PASS>0)&(X.oxoG_PASS>0)&(X.NALT_PASS>0);
tab(passfilt)
XP=trimStruct(X,passfilt)
clf
ALGS={'M1','M2','Strelka1','Strelka2','Snowman','SvABA'}
P=[]
P.log=1
P.angle=25
[h,sv] =  upset_plot(XP,ALGS,P);
set(h(1),'ylim',[0.9 5e6])
AREA='Doc/draft-20201009/plots/'
PSET='REBC_primary_pair_383'
f=[AREA PSET '.Mut_final_algs_UpSet.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')


%%
X1=load('Mutations/REBC_primary_pair_393.maf_ac_pp_ccf_fit_v3.maf.mat')
X1.LABELED_REBC_ID=X1.Tumor_Sample_Barcode;
X1.LABELED_REBC_ID=regexprep(X1.LABELED_REBC_ID,'-YQ-','-');
X1.LABELED_REBC_ID=cellfun(@(x) x(1),regexp(X1.LABELED_REBC_ID,'-TTP','split'));
X1.LABELED_REBC_ID=cellfun(@(x) x(1),regexp(X1.LABELED_REBC_ID,'-TP','split'));

% only analyze pass QC set 383
X1=trimStruct(X1,ismember(X1.LABELED_REBC_ID,T.LABELED_REBC_ID))

tab(X1.LABELED_REBC_ID)
% 
% myhist(X1.t_alt_count,0:50)
% myhist(X.t_alt_count(ismember(X.Variant_Type,{'INS','DEL'})),0:50)  
% tab(X.t_alt_count<1)  % 1248 Strelka 0's 
% printStruct(X,X.t_alt_count<1)


X=trimStruct(X,X.t_ref_count>=0)


X.vaf=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
[i m]=ismember(X.LABELED_REBC_ID,T.LABELED_REBC_ID);

X1.M1=ismember(X1.M1,'True');
X1.M2=ismember(X1.M2,'True');
X1.Strelka1=ismember(X1.Strelka1,'True');
X1.Strelka2=ismember(X1.Strelka2,'True');
X1.Snowman=ismember(X1.Snowman,'True');
X1.SvABA=ismember(X1.SvABA,'True');

figure(3)
clf
ALGS={'M1','M2','Strelka1','Strelka2','Snowman','SvABA'}
P=[]
P.log=1
P.angle=25
[h,sv] =  upset_plot(X1,ALGS,P);
set(h(1),'ylim',[0.9 1e6])
AREA='Doc/draft-20201009/plots/'
PSET='REBC_primary_pair_383'
f=[AREA PSET '.Mut_final_algs_UpSet.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')

figure(4)
[m,k]=sort(-T.WGS_TP_N_nMut)
subplot(4,1,1)
h=bar(T.WGS_TP_N_Case_MEDIAN_COVERAGE(k) ); grid on; ylabel('TP median cov'); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(4,1,2)
h=bar(T.WGS_TP_N_Control_MEDIAN_COVERAGE(k)); grid on; ylabel('N median cov'); ylim([0 150]); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(4,1,3)
h=bar(100*T.WGS_TP_N_purity(k)); grid on; ylabel('Tumor purity (%)'); ylim([0 100]); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(4,1,4)
h=bar(100*T.WGS_TP_N_nMut(k)); grid on; ylabel('Mutations'); set(gca,'yscale','log');ylim([1 5e5]); set(h,'facecolor',0.25*[1 1 1 ])
xlabel('tumors'); set(gca,'ytick',10.^[1:5])
f=[AREA PSET '.Properties.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')


%%
clf
myhist(X1.n_alt_count,0:20,'log')
% conflicting algorithms (eg. M1 and SvABA can't find the same variant)
k=find(X1.M1.*X1.SvABA)
tab(X1.Variant_Type(k))
tab(X1.Genome_Change(k))


k=find(X1.M1);
tab(X1.Variant_Type(k))
k=k(find(ismember(X1.Variant_Type(k),{'INS','DEL'})));

% 140 M1 calls that ended up as indels due to overlap with M2, SvABA, snowman

tab(X1.Variant_Type(k))
tab(X1.M2(k))
tab(X1.Strelka1(k))
tab(X1.Strelka2(k))
tab(X1.Snowman(k))
tab(X1.SvABA(k))


tab(X1.LABELED_REBC_ID(k))
printStruct(X1,k(find(ismember(X1.LABELED_REBC_ID(k),'REBC-ACE3'))))

save(['Mutations/REBC_primary_pair_383.UpSet' today '.mat'],'X','T','X1','PSET')


%%
clear
load('Mutations/REBC_primary_pair_383.UpSet.03Nov2020.mat')
figure(4)
[m,k]=sort(-T.WGS_TP_N_nMut)
subplot(5,1,1)
h=bar(T.WGS_TP_N_Case_MEDIAN_COVERAGE(k) ); grid on; ylabel('TP median cov'); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(5,1,2)
h=bar(T.WGS_TP_N_Control_MEDIAN_COVERAGE(k)); grid on; ylabel('N median cov'); ylim([0 150]); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(5,1,3)
h=bar(100*T.WGS_TP_N_purity(k)); grid on; ylabel('Tumor purity (%)'); ylim([0 100]); set(gca,'xticklabel',[]); set(h,'facecolor',0.25*[1 1 1 ])
subplot(5,1,4)

h=bar(100*T.WGS_TP_N_nMut(k)); grid on; ylabel('Mutations'); set(gca,'yscale','log');ylim([1 5e5]); set(h,'facecolor',0.25*[1 1 1 ])
xlabel('tumors'); set(gca,'ytick',10.^[1:5])
subplot(5,1,5)
h=bar(100*T.WGS_TP_N_nMut(k)); grid on; ylabel('Mutations'); set(gca,'yscale','log');ylim([1 5e5]); set(h,'facecolor',0.25*[1 1 1 ])
xlabel('tumors'); set(gca,'ytick',10.^[1:5])
f=[AREA PSET '.Properties.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')


