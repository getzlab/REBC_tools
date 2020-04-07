function [ CCF, CCFCI,CCFCI95,MODEL,pCCF,multiplicity,NLL,AFM,OVER] = fit_SV_CCF_v3(ALT,REF,NA,NB,NA_CCF,NB_CCF,NCNA,NCNB,purity,PAR)
% fit variant CCF based on 
%   ALT, REF observed allele counts (one at a time) 
%   P tumor purity
%   CN local nromal copy number 
%   CCFC tumor copy number CCF
%   NA, NB minor and major absolute variant copy numbers in tumor
% depth 

if (nargin<10)
    PAR=[];
    PAR.DELTA=0.01;
    PAR.CCF_PRIOR_SLOPE=0.01;
end
CCF_BINS=0:PAR.DELTA:4;
CCF_BINS1=CCF_BINS(CCF_BINS<=1);



CCF_PRIOR=1+0*CCF_BINS; CCF_PRIOR(CCF_BINS>1)=0;  CCF_PRIOR=CCF_PRIOR./sum(CCF_PRIOR);
CCF_PRIOR=(CCF_BINS-0.5)*PAR.CCF_PRIOR_SLOPE*CCF_PRIOR(1)+CCF_PRIOR;
CCF_PRIOR=CCF_PRIOR/mean(CCF_PRIOR(CCF_BINS<=1));CCF_PRIOR(CCF_BINS>1)=0;

% DNA per cell
%tumorDNA=X.SCNA_NA(i).*X.SCNA_NA_CCF(i)+X.SCNA_NB(i).*X.SCNA_NB_CCF(i).*X.SCNA_NCNB(i)+X.SCNA_NCNA(i).*(1-X.SCNA_NA_CCF(i))+X.SCNA_NCNB(i).*(1-X.SCNA_NB_CCF(i));
tumorDNA=NA*NA_CCF+NB*NB_CCF*NCNB+NCNA*(1-NA_CCF)+NCNB*(1-NB_CCF);
CCF=NaN; CCFCI=[NaN NaN];MODEL=0;pCCF=NaN*CCF_BINS1;multiplicity=NaN;NLL=NaN;AFM=NaN;OVER=NaN;CCFCI95=CCFCI;
if tumorDNA<=0
    %no tumor DNA -> no mutation
    fprintf('tumorDNA=0 NA=%d NA_CCF=%.3f NB=%d NB_CCF-%.3f ALT=%d REF=%d',NB,NA_CCF,NB,NB_CCF,ALT,REF)
    return;
end
normalDNA=NCNA+NCNB;
totalDNA = tumorDNA*purity + normalDNA*(1-purity);

X.dna_fraction_in_tumor=(purity*tumorDNA)./(purity*tumorDNA+(1-purity)*normalDNA);

% skip events that ABSOLUTE skipped
if isnan(X.dna_fraction_in_tumor), return; end
%
% six models: mutation on A or B, on SCNA, normal, or both
% ALT DNA per cell
ALT_DNA=zeros(6,length(CCF_BINS));
ALT_LIKE_FULL=zeros(6,length(CCF_BINS));
ALT_LIKE_WEIGHT=zeros(6,length(CCF_BINS));
%ALT_LIKE=zeros(6,length(CCF_BINS));
% limit ccf to SCNA CCF of whatever alle mutaiton is on
% 1) mutation on A allele NA, m=1
ccfmax=1/PAR.DELTA+1+zeros(6,1);
ccf1=CCF_BINS; ccfmax(1)=NA_CCF;
ALT_DNA(1,:) = ccf1*purity*(NA_CCF)*(NA>0);
% 2) mutation on B allele NB, m=1
ccf1=CCF_BINS; ccfmax(2)=NB_CCF;
ALT_DNA(2,:) = ccf1*purity*(NB_CCF)*(NB>0);
% 3) mutation on A allele NCNA, m=1
ccf1=CCF_BINS; ccfmax(3)=1-NA_CCF;
ALT_DNA(3,:) = ccf1*purity*(NCNA*(1-NA_CCF));
% 4) mutation on B allele NCNB, m=1
ccf1=CCF_BINS; ccfmax(4)=1-NB_CCF;
ALT_DNA(4,:) = ccf1*purity*(NCNB*(1-NB_CCF));
% 5) mutation on A allele NCNA and NA m=NA
ccf1=CCF_BINS; ccfmax(5)=1;
if (NA>0)
    ALT_DNA(5,:) = ccf1*purity*( (NA*NA_CCF) + (NCNA*(1-NA_CCF)) );
end
% 6) mutation on B allele NCNB and NB m=NB
ccf1=CCF_BINS; ccfmax(6)=1;
if (NB>0)
    ALT_DNA(6,:) = ccf1*purity*( (NB*NB_CCF) + (NCNB*(1-NB_CCF)) );
end
like1=0;
model1=0;
af1=0;
P1=0*ccf1;
for m1=1:6
    AF_DNA(m1,:)=ALT_DNA(m1,:)./totalDNA;
    AF_DNA(m1,AF_DNA(m1,:)>1)=0;
    % likelihood for ccfs + extended ccfs to allow for binominal AF fluctuation above max CCF
    ALT_LIKE_FULL(m1,:)=binopdf(ALT,ALT+REF,AF_DNA(m1,:));
    % stdev
    [mu1,sigma1]=stats_mu_sigma(CCF_BINS,ALT_LIKE_FULL(m1,:) );
    % weight
    W=exp(-0.5*(CCF_BINS-ccfmax(m1,:)).^2./sigma1^2);
    W=W/max(W);
    W(CCF_BINS<=ccfmax(m1))=1;
    % restrict to possible ccf range
    % prob CCF conditional on CCF < ccfmax - attenutate ALT_LIKE
    P=ALT_LIKE_FULL(m1,:).*W;
    PX=P(ccf1>ccfmax(m1));
    kend=find(ccf1<=ccfmax(m1)); kend=kend(end);
    if (sum(PX)>0)
        P(kend)=P(kend)+sum(PX);
        P((kend+1):end)=0;
    end
    P=P.*CCF_PRIOR;
    ALT_LIKE_WEIGHT(m1,:)=P;
    % pick top likelihood - possible range
    if max(P)>like1
        [like1,k1]=max(P);
        model1=m1;
        af1=AF_DNA(m1,k1);
        ccfmax1=ccfmax(m1);
        P1=P;
    end
end
% possible ccf probabilties - without weight W
pccf=ALT_LIKE_FULL(model1,:);  pccf=pccf/sum(pccf);
pccfx=sum(pccf(ccf1>ccfmax1));
kend=find(ccf1<=ccfmax1); kend=kend(end);
if (pccfx>0)
    pccf(kend)=pccf(kend)+pccfx;
    pccf((kend+1):end)=0;
end

if isnan(sum(pccf))
   fprintf('pccf NaN ALT=%d, REF=%d, P',ALT,REF,P,CN,NA,NB,NA_CCF,NB_CCF,purity)
   return;
end

NLL=-log10(like1);
MODEL=model1;
AFM=af1;
OVER=pccfx/sum(pccf);



pCCF=pccf(CCF_BINS<=1)/sum(pccf(CCF_BINS<=1));

%X.ccf_hat(i,1)=sum(ccf.*X.pccf(i,:))./sum(X.pccf(i,:));

nb=length(CCF_BINS1);
cccf=cumsum(pCCF);
ilow=sum(cccf<0.16);
ihigh=sum(cccf<0.84)+1;
if (ilow<1), ilow=1;end
if (ihigh<=ilow), ihigh=ilow+1;end
if (ihigh>nb), ihigh=nb;end
ccf_CI_low=CCF_BINS(ilow);
ccf_CI_high=CCF_BINS(ihigh);

ilow=sum(cccf<0.025);
ihigh=sum(cccf<0.975)+1;
imedian=sum(cccf<0.5);
if (imedian<1), imedian=1; end
[~,imax]=max(pCCF);
if (ilow<1), ilow=1;end
if (ihigh<=ilow), ihigh=ilow+1;end
if (ihigh>nb), ihigh=nb;end
ccf_CI95_low=CCF_BINS(ilow);
ccf_CI95_high=CCF_BINS(ihigh);
ccf_mode=CCF_BINS(imax(1));
ccf_mean=sum(CCF_BINS1.*pCCF)./sum(pCCF);
ccf_median=CCF_BINS(imedian(1));
multiplicity=1;
if model1==5
    multiplicity=NA;
end
if model1==6
    multiplicity=NB;
end

CCF =ccf_median;
CCFCI=[ccf_CI_low ccf_CI_high];
CCFCI95=[ccf_CI95_low ccf_CI95_high];

end

function [mu,sig]=stats_mu_sigma(x,p)
 p=p/sum(p);
 mu=sum(x.*p);
 sig=sqrt( sum(x.^2.*p) - mu.^2);
end

%     
%     
% DEP=ALT+REF;
% % vary CCFV from 0 to 1 in NB bins
% nb=300;DCCF=1/nb;
% CCFV=DCCF:DCCF:1;
% % tumor alt absolute copy number 
% CT=NA+NB;
% 
% % anomyous matlab functions for ALT and REF counts 
% 
% ALT1= @(ccf) DEP*P*ccf;
% ALT21= @(ccf) DEP*P*ccf;
% ALT2a= @(ccf) DEP*P*((ccf-CCFC) + CCFC*NA);
% ALT2b= @(ccf) DEP*P*((ccf-CCFC) + CCFC*NB);
% ALT3= @(ccf) P*ccf;
% 
% REF1=  @(ccf) DEP* ( P* (ccf*(CN-1)+(1-ccf-CCFC)*CN + CCFC*CT) + (1-P)*CN) ;
% REF21= @(ccf) DEP* ( P* (ccf*(1-CCFC)*(CN-1)+(1-ccf)*CN + CCFC*(CT-1)) + (1-P)*CN) ;
% REF2a= @(ccf) DEP* ( P* (ccf*(1-CCFC)*(CN-1)+(1-ccf)*CN + CCFC*(CT-NA)) + (1-P)*CN) ;
% REF2b= @(ccf) DEP* ( P* (ccf*(1-CCFC)*(CN-1)+(1-ccf)*CN + CCFC*(CT-NB)) + (1-P)*CN) ;
% %REF3=  @(ccf) DEP* ( P* ((CCFC*CT-ccf)+(1-CCFC)*CN) + (1-P)*CN) ;
% %REF3=  @(ccf) DEP* ( P* ( (1-CCFC)*CN - ccf + CCFC*(1-ccf)*CT ) + (1-P)*CN) ;
% REF3=  @(ccf) DEP* ( P* ( (1-CCFC)*CN - CCFC*(1-ccf)*CT  + ccf*(CT-1) ) + (1-P)*CN) ;
% 
% alt=zeros(size(CCFV));
% ref=zeros(size(CCFV));
% 
% for i=1:nb
%     alt1(i)  = ALT1(CCFV(i));
%     alt21(i) = ALT21(CCFV(i));
%     alt2a(i) = ALT2a(CCFV(i));
%     alt2b(i) = ALT2b(CCFV(i));
%     alt3(i)  = ALT3(CCFV(i));
%     % 
%     ref1(i)  = REF1(CCFV(i));
%     ref21(i) = REF21(CCFV(i));
%     ref2a(i) = REF2a(CCFV(i));
%     ref2b(i) = REF2b(CCFV(i));
%     ref3(i)  = REF3(CCFV(i));
% 
% end
% EPS=1e-200; % eps(1e-150)
% palt1 =  EPS+poisspdf(ALT,alt1);
% palt21 = EPS+poisspdf(ALT,alt21);
% palt2a = EPS+poisspdf(ALT,alt2a);
% palt2b = EPS+poisspdf(ALT,alt2b);
% palt3 =  EPS+poisspdf(ALT,alt3);
% 
% pref1 =  EPS+poisspdf(REF,ref1);
% pref21 = EPS+poisspdf(REF,ref21);
% pref2a = EPS+poisspdf(REF,ref2a);
% pref2b = EPS+poisspdf(REF,ref2b);
% pref3 =  EPS+poisspdf(REF,ref3);
% 
% p1=pref1.*palt1;
% p21=pref21.*palt21;
% p2a=pref2a.*palt2a;
% p2b=pref2b.*palt2b;
% p3=pref3.*palt3;
% 
% p1(isnan(p1))=0;
% p2a(isnan(p2a))=0;
% p2b(isnan(p2b))=0;
% p3(isnan(p3))=0;
% 
% % variant can't be on deleted DNA
% if (CT==0)
%     p21=0*p1;
%     p2a=0*p2a;
%     p2b=0*p2b;
% end  
% 
% [px(1) k(1)]=max(p1);
% [px(2) k(2)]=max(p21);
% [px(3) k(3)]=max(p2a);
% [px(4) k(4)]=max(p2b);
% [px(5) k(5)]=max(p3);
% 
% if all(isnan(p1)); px(1)=0;    end
% if all(isnan(p21)); px(2)=0;    end
% if all(isnan(p2a)); px(3)=0;    end
% if all(isnan(p2b)); px(4)=0;    end
% if all(isnan(p3)); px(5)=0;    end
% 
% [PX,K]=max(px);
% 
% k1=k(K);
% CCF=CCFV(k1);
% MODEL = K;
% pCCF=p1; m=1;
% if (K==2), pCCF=p21; m=1; end
% if (K==3), pCCF=p2a; m=NA; end
% if (K==4), pCCF=p2b; m=NB; end
% if (K==5), pCCF=p3; m=1; end
% 
% pCCF=pCCF/sum(pCCF);
% cCCF=cumsum(pCCF);
% klow=find(cCCF>=0.05); klow=klow(1);
% khigh=find(cCCF>=0.95); khigh=khigh(end);
% if (klow<1), klow=1; end
% if (khigh>nb), khigh=nb; end
% CCFCI(1)=CCFV(klow);
% CCFCI(2)=CCFV(khigh);
% NLL=-log10(PX);
% 
% end
%%
function test1

NA=0
NB=2
NA_CCF=1
NB_CCF=0.95
ALT=20
REF=50
purity=0.6
NCNA=1
NCNB=1

 [ CCF, CCFCI,CCFCI95,MODEL,pCCF,multiplicity,NLL,AFM,OVER] = fit_SV_CCF_v3(ALT,REF,NA,NB,NA_CCF,NB_CCF,NCNA,NCNB,purity)
end

