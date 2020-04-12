function X=SV_CCF_v3(X1,ST1,gender,PAR)
% function X=SV_CCF_V3(X1,ST1,WINDOW)
%
%
if (nargin<4)
    PAR=[];
    PAR.DELTA=0.01;
    PAR.CCF_PRIOR_SLOPE=0.01;
    PAR.ASSUME_NORMAL_X_WHEN_MISSING=1;
    PAR.WINDOW=10e3;
    PAR.CLONAL_CCF_THRESHOLD=0.75;
end
CCF_BINS=0:PAR.DELTA:1;
NB=length(CCF_BINS);

N=length(X1.chr1);
X1.x1=xhg19(X1.chr1,X1.pos1);
X1.x2=xhg19(X1.chr2,X1.pos2);
X1.VCF_VAF= NaN*X1.pos1;
X1.VCF_VAF_95CI_low = NaN*X1.pos1;
X1.VCF_VAF_95CI_high = NaN*X1.pos1;
X1.VCF_VAF_CI_low = NaN*X1.pos1;
X1.VCF_VAF_CI_high = NaN*X1.pos1;

X1.purity=NaN*X1.pos1;
X1.ploidy=NaN*X1.pos1;

X1.SCNA_NA1=NaN*X1.pos1;
X1.SCNA_NB1=NaN*X1.pos1;
X1.SCNA_NA_CCF1=NaN*X1.pos1;
X1.SCNA_NB_CCF1=NaN*X1.pos1;

X1.SCNA_NA2=NaN*X1.pos1;
X1.SCNA_NB2=NaN*X1.pos1;
X1.SCNA_NA_CCF2=NaN*X1.pos1;
X1.SCNA_NB_CCF2=NaN*X1.pos1;

X1.ccf_hat=NaN*X1.pos1;
X1.ccfCI_low=NaN*X1.pos1;
X1.ccfCI_high=NaN*X1.pos1;
X1.ccfMultiplicity=NaN*X1.pos1;
X1.ccfNegativeLogLikelihood=NaN*X1.pos1;
    
X1.ccf_hat_break1=NaN*X1.pos1;
X1.ccfMultiplicity_break1=NaN*X1.pos1;
X1.ccfNegativeLogLikelihood_break1=NaN*X1.pos1;
X1.ccf_hat_break2=NaN*X1.pos1;
X1.ccfMultiplicity_break2=NaN*X1.pos1;
X1.ccfNegativeLogLikelihood_break2=NaN*X1.pos1;

X1.ccf_hat_break1=NaN*X1.pos1;
X1.ccfMultiplicity_break1=NaN*X1.pos1;
X1.ccfNegativeLogLikelihood_break1=NaN*X1.pos1;
X1.ccf_hat_break2=NaN*X1.pos1;
X1.ccfMultiplicity_break2=NaN*X1.pos1;
X1.ccfNegativeLogLikelihood_break2=NaN*X1.pos1;

X1.clonal=NaN*X1.pos1;

if N<1
    X=X1;
    return
end

X1.VCF_VAF=2*X1.VCF_TALT./(X1.VCF_TALT+X1.VCF_TREF);
[X1.VCF_VAF, ci]=binofit(X1.VCF_TALT*2,(2*X1.VCF_TALT+X1.VCF_TREF),0.05);
X1.VCF_VAF_95CI_low = ci(:,1);
X1.VCF_VAF_95CI_high = ci(:,2);
[~, ci]=binofit(2*X1.VCF_TALT,(2*X1.VCF_TALT+X1.VCF_TREF),0.32);
X1.VCF_VAF_CI_low = ci(:,1);
X1.VCF_VAF_CI_high = ci(:,2);
    
ST1.x1=xhg19(ST1.Chromosome,ST1.Start_bp);
ST1.x2=xhg19(ST1.Chromosome,ST1.End_bp);

ST1.CNT=(ST1.tau-ST1.CCN.*(1-ST1.purity))./ST1.purity;
ST1.AF2CCF=(ST1.CCN.*(1-ST1.purity) + ST1.purity.*ST1.CNT)./ST1.purity;

purity=ST1.purity(1);
ploidy=ST1.ploidy(1);

for i1=1:N
    k1b=find((X1.x1(i1)>=ST1.x1)&(X1.x1(i1)<=ST1.x2));
    k2b=find((X1.x2(i1)>=ST1.x1)&(X1.x2(i1)<=ST1.x2));

    % which side of break to check aSCNAs 
    % forward: str=0, should check left D = -1
    % reverse: str=1, should check right D= +1
    % wrong way 07Apr2020: D1=PAR.WINDOW*(1-2*X1.str1(i1)); D2=PAR.WINDOW*(1-2*X1.str2(i1)); 
    D1=PAR.WINDOW*(2*X1.str1(i1)-1);
    D2=PAR.WINDOW*(2*X1.str2(i1)-1);
 
    x1a=X1.x1(i1)+D1;
    x2a=X1.x2(i1)+D2;
    k1=find((x1a>=ST1.x1)&(x1a<=ST1.x2));
    k2=find((x2a>=ST1.x1)&(x2a<=ST1.x2));
    if (k1~=k1b)
        fprintf('%d\t%d\n',[ST1.CNT(k1) ST1.CNT(k1b)])
        fprintf('%d\t%d\n',[ST1.tau(k1) ST1.tau(k1b)])
        printStruct(X1,i1)
    end
    if (k2~=k2b)
        fprintf('%d\t%d\n',[ST1.CNT(k2) ST1.CNT(k2b)])
        fprintf('%d\t%d\n',[ST1.tau(k2) ST1.tau(k2b)])
        printStruct(X1,i1)
    end
    if (~isempty(k1))&&(~isnan(ST1.NA_CCF(k1)))
        X1.SCNA_NA1(i1,1)=ST1.NA(k1);
        X1.SCNA_NB1(i1,1)=ST1.NB(k1);
        X1.SCNA_NA_CCF1(i1,1)=ST1.NA_CCF(k1);
        X1.SCNA_NB_CCF1(i1,1)=ST1.NB_CCF(k1);
    end
    if (~isempty(k2))&&(~isnan(ST1.NA_CCF(k2)))
        X1.SCNA_NA2(i1,1)=ST1.NA(k2);
        X1.SCNA_NB2(i1,1)=ST1.NB(k2);
        X1.SCNA_NA_CCF2(i1,1)=ST1.NA_CCF(k2);
        X1.SCNA_NB_CCF2(i1,1)=ST1.NB_CCF(k2);
    end
end

ccf_hat(2)=NaN; CCFCI(2,2)=NaN;CCFCI95(2,2)=NaN;MODEL(2)=NaN;pCCF(2,NB)=NaN;multiplicity(2)=NaN;NLL(2)=NaN;AFM(2)=NaN;OVER(2)=NaN;

for i=1:N
    
    X1.purity(i)=purity;
    X1.ploidy(i)=ploidy;
    NCNA1=1; NCNB1=1;
    NCNA2=1; NCNB2=1;
    if PAR.ASSUME_NORMAL_X_WHEN_MISSING&&ismember(upper(gender),{'M','MALE'})
       if(X1.chr1==23),NCNA1=0; end
       if(X1.chr2==23),NCNA2=0; end
    end            
    % which end for SCNA info NCN CCF_SCAN NA NB ? 
    %[ CCF, CCFCI,CCFCI95,MODEL,pCCF,multiplicity,NLL,AFM,OVER] = fit_SV_CCF_v3(ALT,REF,NA,NB,NA_CCF,NB_CCF,NCNA,NCNB,purity)
    [ ccf_hat(1), CCFCI(1,:),CCFCI95(1,:),MODEL(1),pCCF(1,:),multiplicity(1),NLL(1),AFM(1),OVER(1)] = fit_SV_CCF_v3(2*X1.VCF_TALT(i),X1.VCF_TREF(i),X1.SCNA_NA1(i),X1.SCNA_NB1(i),X1.SCNA_NA_CCF1(i),X1.SCNA_NB_CCF1(i),NCNA1,NCNB1,purity);
    [ ccf_hat(2), CCFCI(2,:),CCFCI95(2,:),MODEL(2),pCCF(2,:),multiplicity(2),NLL(2),AFM(2),OVER(2)] = fit_SV_CCF_v3(2*X1.VCF_TALT(i),X1.VCF_TREF(i),X1.SCNA_NA2(i),X1.SCNA_NB2(i),X1.SCNA_NA_CCF2(i),X1.SCNA_NB_CCF2(i),NCNA2,NCNB2,purity);
    if isnan(ccf_hat(1))
        pccf=pCCF(2,:);
    end
    if isnan(ccf_hat(2))
        pccf=pCCF(1,:);
    end    
    if prod(ccf_hat)>0
        pccf=(pCCF(1,:)+pCCF(2,:))/2;
    end
    if isnan(ccf_hat(1))&isnan(ccf_hat(2))
        continue;
    end    
    cccf=cumsum(pccf);
    imedian=sum(cccf<0.5);
    if (imedian<1), imedian=1; end
    [~,~]=max(pccf);
    X1.ccf_hat(i)=CCF_BINS(imedian(1));
    ilow=sum(cccf<0.16);
    ihigh=sum(cccf<0.84)+1;
    if (ilow<1), ilow=1;end
    if (ihigh<=ilow), ihigh=ilow+1;end
    if (ihigh>NB), ihigh=NB;end
    X1.ccfCI_low(i)=CCF_BINS(ilow);
    X1.ccfCI_high(i)=CCF_BINS(ihigh);
    
    X1.ccfMultiplicity(i)=nanmin(multiplicity);
    X1.ccfNegativeLogLikelihood(i)=nanmax(NLL);
    
    X1.ccf_hat_break1(i)=ccf_hat(1);
    X1.ccfMultiplicity_break1(i)=multiplicity(1);
    X1.ccfNegativeLogLikelihood_break1(i)=NLL(1);
    X1.ccf_hat_break2(i)=ccf_hat(2);
    X1.ccfMultiplicity_break2(i)=multiplicity(2);
    X1.ccfNegativeLogLikelihood_break2(i)=NLL(2);
    
    X1.clonal(i)=1*(X1.ccf_hat(i)>=PAR.CLONAL_CCF_THRESHOLD);
end
X=X1;
end

%%
function test
%%
clear
close all
addpath('/Volumes/GoogleDrive/My Drive/Cancer/tools/matlab')
addpath('/Volumes/GoogleDrive/My Drive/Cancer/REBC/tools/matlab')

PSET='REBC_primary_pair_393'; 
%PSET='REBC-NT-NB_240'; 
%DAY='17Dec2019'
DAY='20Feb2020'

AREA=['/Volumes/GoogleDrive/My Drive/Cancer/REBC/SV/']
cd(AREA)
X=load_tsv([AREA PSET '.BP.vcf.filter.NALT.TALT.' DAY '.tsv'])

ST=load_tsv(['../CN/CCF/' PSET '.alleliccapseg_pp_ccf_fit_v3_collapsed_seg.tsv'])


P=load_tsv('../FC/REBC.pairs.27Nov2019.tsv')
PS=load_tsv(['../FC/' PSET '.pair_set.tsv'])
P=RenameField(P,'entity_pair_id','pair_id')
P=trimStruct(P,ismember(P.pair_id,PS.pair))
PG=load_tsv('../FC/load_REBC_sample_gender.tsv')
[i,m]=ismember(P.case_sample,PG.update_sample_id)
tab(i)
P.Gender(i,1)=PG.Gender(m(i))


% only calls in P
[i m]=ismember(X.individual, P.pair_id);
tab(i)
X=trimStruct(X,i)
[i m]=ismember(X.individual, P.pair_id)
tab(i)

pid=unique(X.individual)
NP=length(pid)
[i m]=ismember(ST.sample,P.pair_id);
tab(i)

x1=[]
V=[]
for i=1:NP
    pid1=pid{i};
    P1=trimStruct(P,ismember(P.pair_id,pid{i}));
    X1=trimStruct(X,ismember(X.individual,pid1));
    ST1=trimStruct(ST,ismember(ST.sample,pid1));
    fprintf('%s, %d\n',pid1,length(X1.chr1));
    
    x1=SV_CCF_v3(X1,ST1,P1.Gender);
    
    if isempty(V)
        V=x1;
    else
        V=mergeStruct(V,x1);
    end
end


PAR=[]

PAR.YLAB='Fit CCF'
PAR.XLAB='VAF'
hist2(V.VCF_VAF,V.ccf_hat,50,50,[],PAR)
    
printStruct(V,-1,[AREA PSET '.BP.vcf.filter.NALT.TALT.CCF.V3.' TODAY '.txt'])
save([AREA PSET '.BP.vcf.filter.NALT.TALT.CCF.V3.' TODAY '.mat'],'-struct','V')
end
%%

function test1
%%
clear
PSET='REBC_primary_pair_393'; 
DAY='20Feb2020'
AREA=['/Volumes/GoogleDrive/My Drive/Cancer/REBC/SV/']
V=load([AREA PSET '.BP.vcf.filter.NALT.TALT.CCF.V3.' TODAY '.mat'])
k=find((V.ccf_hat<0.05)&(V.VCF_VAF>0.1))
tab(~isnan(V.ccf_hat))

[n,b]=histc(V.purity,0:0.1:1)
for p=1:10
    k=find(b==p);
    plot(V.VCF_VAF(k),V.ccf_hat(k),'o')
    hold on
end
hold off
legend(regexp(num2str((1:10)/10),' +','split'))

k=find((V.ccf_hat<(V.VCF_VAF*1.25)))
hold on
plot(V.VCF_VAF(k),V.ccf_hat(k),'kx')
hold off

end

function make1
main='MAF_AC_PP_CCF_fit_v3'
flist = matlab.codetools.requiredFilesAndProducts( main );

mkdir(['/tmp/' main])
cd(['/tmp/' main])
for i=1:length(flist)
    cmd=['scp  "' flist{i} '" .']
    unix(cmd)
end
unix('ls -lath')
end

