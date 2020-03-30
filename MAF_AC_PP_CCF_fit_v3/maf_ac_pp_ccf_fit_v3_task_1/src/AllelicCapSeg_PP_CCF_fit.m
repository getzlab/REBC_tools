function [X]=AllelicCapSeg_PP_CCF_fit(AC,purity,ploidy,gender,PAR)
%[X]=AllelicCapSeg_PP_CCF_fit (AC,purity,ploidy,PAR)
%
%  AllelicCapSeg_PP_CCF_fit adds estimated tumor SCNA allele copy number [NA,NB] and 
%  CCF's to AC based on puity and ploidy (from ABSOLUTE or other method).
%  This function makes estimates for all segments, even X and Y when gender
%  info is available from PAR. 
%  AllelicCapSeg_PP_CCF_fit operates on one sample at a time (one purity, ploidy, gender) 
%
%  Note: X output roughly correpsonds to the ABSOLUTE SEG_TAB file, with a few caveats
%      1) The model includes 1 SCNA 'state' (NA,NB) and 1 normal state (1,1) 
%      2) The one SCNA state has one CCF and the normal CN occurs with fraction (1-CCF) 
%      3) includes X and Y Chromosome estimates
%      4) includes [NA NB] and CCF even for segments w/o het (mu_major/minor) info 
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
%       CN_CHISQ: Bias toward lower CN solutions added to chisquare (0.5)
%       HIGH_CCF_CHISQ: Bias toward higher CCF solutions added to chisquare (1)
%       SMOOTH_SPAN: taps in chisquare smoothing function (7)
%       CCF_BINS: bins over CCF range (0-1)
%       SIGMA_THRESHOLD: signal:noise threshold for SCNA (2)
%
%  output: 
%     X: AC struct with new fields: 
%       NA: SCNA Allele A copy number (lower or equal CN)
%       NB: SCNA Allele B copy number (higher CN)
%       chisq: fit chisquared
%       NA_nextbest: next best NA
%       NB_nextbest: next best NB
%       chisq_nextbest: next best chisq
%       CCF_prob: vector of CCF probalities
%       CCF_hat: CCF estimate (mode)
%       CCF_Low: CCF 95% low CI
%       CCF_High: CCF 95% high CI
%       IS_SCNA: tumor SCNA flag
%       CCF_prob_nextbest: next best vector of CCF probalities
%       CCF_hat_nextbest: next best CCF estimate (mode)
%       purity: 
%       ploidy: 
%       gender: MALE or FEMALE
%       CNN: normal CN state (2, except 1 for male X, Y, 0 for female Y)
%       muA_fit: fit value of mu_minor
%       muB_fit: fit value of mu_major
%       tau_fit: fit value of tau
%       muA_nextbest: second best fit value of mu_minor
%       muB_nextbest: second best fit value of mu_major
%       tau_nextbest: second best fit value of tau
%
% Amaro / Chip 11 Mar 2015
%
if (nargin<5)
    PAR=[];
    PAR.CCF_BINS=0:0.01:1;
    PAR.SMOOTH_SPAN=7;
    PAR.HIGH_CCF_CHISQ=0.5;
    PAR.CN_CHISQ=0.5; 
    PAR.MAXCN=15;
    PAR.SIGMA_MIN=0.05;
    PAR.SIGMA_THRESHOLD=2;
    PAR.X_MALE_TAU_SCALE = 1
end
if (nargin<4)
    msgID = 'AllelicCapSeg_PP_CCF_fit:Input';
    msg = 'Not enough inputs: AC, purity, ploidy, gender';
    throw(MException(msgID,msg));
end
if ~isfield(PAR,'CCF_BINS')
    PAR.CCF_BINS=0:0.01:1; 
end
if ~isfield(PAR,'SMOOTH_SPAN')
    PAR.SMOOTH_SPAN=7; 
end
if ~isfield(PAR,'HIGH_CCF_CHISQ')
    PAR.HIGH_CCF_CHISQ=1; 
end
if ~isfield(PAR,'CN_CHISQ')
    PAR.CN_CHISQ=0.5; 
end
if ~isfield(PAR,'MAXCN')
    PAR.MAXCN=15; 
end
if ~isfield(PAR,'SIGMA_MIN')
    PAR.SIGMA_MIN=0.05; 
end
if ~isfield(PAR,'SIGMA_THRESHOLD')
    PAR.SIGMA_THRESHOLD=2;
end

CCF=PAR.CCF_BINS;

% set nan sigmas to something other than nan
k=find(isnan(AC.sigma_major))
AC.sigma_minor(k)=nanmedian(AC.sigma_minor)
AC.sigma_major(k)=nanmedian(AC.sigma_major)
% intitialize AC fields
AC.sigma_major=sqrt(AC.sigma_major.^2+PAR.SIGMA_MIN.^2);
AC.sigma_minor=sqrt(AC.sigma_minor.^2+PAR.SIGMA_MIN.^2);
AC.sigma_tau=sqrt(AC.sigma_tau.^2+PAR.SIGMA_MIN.^2);
AC.NA=NaN*AC.length;
AC.NB=NaN*AC.length;
AC.chisq=NaN*AC.length;
AC.CCF_prob=NaN(length(AC.length),length(CCF));
AC.CCF_hat=NaN*AC.length;
AC.CCF_Low=NaN*AC.length;
AC.CCF_High=NaN*AC.length;
AC.muA_fit=NaN*AC.length;
AC.muB_fit=NaN*AC.length;
AC.tau_fit=NaN*AC.length;
AC.NA_nextbest=NaN*AC.length;
AC.NB_nextbest=NaN*AC.length;
AC.chisq_nextbest=NaN*AC.length;
AC.CCF_prob_nextbest=NaN(length(AC.length),length(CCF));
AC.CCF_hat_nextbest=NaN*AC.length;
AC.muA_nextbest=NaN*AC.length;
AC.muB_nextbest=NaN*AC.length;
AC.tau_nextbest=NaN*AC.length;
AC.purity=NaN*AC.length;
AC.ploidy=NaN*AC.length;
AC.gender=repmat(gender,length(AC.length),1);


AC.IS_SCNA=NaN*AC.length;
XCN=2; % X normal CN 
if ismember(upper(gender),{'M','MALE'})
    XCN=1;
    
    k=find(AC.Chromosome>=23)
    AC.tau(k)=AC.tau(k)/PAR.X_MALE_TAU_SCALE;
    AC.mu_minor(k)=0;
    AC.mu_major(k)=AC.tau(k)
    
end


% loop over SCNA segments
for i=1:length(AC.sigma_major)
    NB=length(PAR.CCF_BINS);
    chi=repmat(inf,[PAR.MAXCN+1 PAR.MAXCN+1 NB]);
    chisq=repmat(inf,[PAR.MAXCN+1 PAR.MAXCN+1]);
    chi1=inf;
    NA1=NaN;
    NB1=NaN;
    tau1=NaN;
    muA1=NaN;
    muB1=NaN;
    chi2=NaN;
    NA2=NaN;
    NB2=NaN;
    tau2=NaN;
    muA2=NaN;
    muB2=NaN;
    chix1=0;
    chix2=0;
    NCN=2;
    if (AC.Chromosome(i)==23)
        NCN=XCN;
    end
    if (AC.Chromosome(i)==24)
        NCN=1-XCN;
    end    
    if (NCN<1), continue; end  % no Y

    for NA=0:PAR.MAXCN
        for NB=NA:PAR.MAXCN
            if (NA==1)&(NB==1)&(NCN==2), continue; end  % not an SCNA 
            if (NA==0)&(NB==1)&(NCN==1), continue; end  % not an SCNA 
            if isnan(AC.mu_major(i)),    continue; end  % no way to determine NA NB locally 
            muA=(purity*CCF*(NA-1)+1)/(purity*(ploidy/2)+(1-purity));
            muB=(purity*CCF*(NB-1)+1)/(purity*(ploidy/2)+(1-purity));
            %muA=((purity*(CCF*(NA-1)+1))+1-purity)/(purity*(ploidy/2)+(1-purity));
            %muB=((purity*(CCF*(NB-1)+1))+1-purity)/(purity*(ploidy/2)+(1-purity));
            %muA=purity*(CCF*(NA-1)+1)+1-purity;
            %muB=purity*(CCF*(NB-1)+1)+1-purity;
            if NCN==2
                tau=2*(purity*CCF*(NA+NB-2)+2)/(purity*ploidy+2*(1-purity));
                %tau=2*((purity*(CCF*(NA+NB-2)+2)+2*(1-purity))/(purity*ploidy+2*(1-purity)));
                %tau=purity*(CCF*(NA+NB-2)+2)+2*(1-purity);
                chi(NA+1,NB+1,:)=((AC.mu_minor(i)-muA)/AC.sigma_minor(i)).^2+((AC.mu_major(i)-muB)/AC.sigma_major(i)).^2;
            else
                if NA>0, continue; end
                tau=muB;
                %tau=NCN*((purity*(CCF*(NA+NB-NCN)+NCN)+NCN*(1-purity))/(purity*ploidy+NCN*(1-purity)));
                chi(NA+1,NB+1,:)=((AC.tau(i)-tau)/AC.sigma_tau(i)).^2;
            end
            % single allele SCNA bias when not doubled           
            % dont add junk to chi
            %if (NB>0)& ( ((NA==NB)&(ploidy<3)) |  ( (NA~=1)&(NB~=1)) )
            %    chi(NA+1,NB+1,:)= chi(NA+1,NB+1,:)  + (PAR.CN_CHISQ+((AC.tau(i)-NCN)/AC.sigma_tau(i))^-2);
            %end
            % high CCF bias = low CCF penalty
            chi(NA+1,NB+1,:)= squeeze(chi(NA+1,NB+1,:))'  + (1-CCF)*PAR.HIGH_CCF_CHISQ ;
            % far from normal NA+NB penalty
            chi(NA+1,NB+1,:)= squeeze(chi(NA+1,NB+1,:))'  + abs(NA+NB-2)*PAR.CN_CHISQ ;
            %smooth to lessen impact of coincidentially sharp local minima due to CCF binning
            chi(NA+1,NB+1,:)=smooth(squeeze(chi(NA+1,NB+1,:)),PAR.SMOOTH_SPAN);
            chisq(NA+1,NB+1)=min(chi(NA+1,NB+1,:),[],3);
            
            if chisq(NA+1,NB+1)<chi1
                % stash second best solution
                NA2=NA1;
                NB2=NB1;
                chi2=chi1;
                muA2=muA1;
                muB2=muB1;
                tau2=tau1;
                % best solution
                NA1=NA;
                NB1=NB;
                [chi1,k1]=min(chi(NA+1,NB+1,:),[],3);
                muA1=muA(k1);
                muB1=muB(k1);
                tau1=tau(k1);
            elseif chisq(NA+1,NB+1)<chi2
                NA2=NA;
                NB2=NB;
                [chi2,k2]=min(chi(NA+1,NB+1,:),[],3);
                muA2=muA(k2);
                muB2=muB(k2);
                tau2=tau(k2);
            end
            
        end
    end
    
    % Amaros change of mind about best and next best to favor NA~=NB solutions
    if (NA1==NB1)&(NA2~=NB2)&(NB1>0)&((chi2-chi1)<PAR.CN_CHISQ)&(ploidy<3)
        NA1T=NA1;        NB1T=NB1;        chi1T=chi1;        
        muA1T=muA1;      muB1T=muB1;      tau1T=tau1;        
        NA1=NA2;         NB1=NB2;          chi1=chi2;        
        muA1=muA2;       muB1=muB2;        tau1=tau2;        
        NA2=NA1T;        NB2=NB1T;        chi2=chi1T;        
        muA2=muA1T;      muB2=muB1T;      tau2=tau1T;        
    end
    
    if ~isnan(NA1+NB1)
       L=squeeze(exp(-chi(NA1+1,NB1+1,:)/2));
       CL=cumsum(L)/sum(L);
    else
        L=NaN;
        CL=NaN;
    end
    if ~isnan(NA2+NB2)
       Lnext=squeeze(exp(-chi(NA2+1,NB2+1,:)/2));
       CLnext=cumsum(Lnext)./sum(Lnext);
    else
        Lnext=NaN;
        CLnext=NaN;
    end
    
    AC.NA(i,1)=NA1;
    AC.NB(i,1)=NB1;
    AC.chisq(i,1)=chi1;
    AC.NA_nextbest(i,1)=NA2;
    AC.NB_nextbest(i,1)=NB2;
    AC.chisq_nextbest(i,1)=chi2;
    AC.CCF_prob(i,:)=squeeze(L)./sum(L(:));
    if ~isnan(NA1+NB1)  
        AC.CCF_hat(i,1)=CCF(find(squeeze(chi(NA1+1,NB1+1,:))==min(squeeze(chi(NA1+1,NB1+1,:))),1,'first'));
    end
    AC.CCF_Low(i,1)=CCF(max([find(CL<=.025,1,'Last') 1]));
    AC.CCF_High(i,1)=CCF(min([find(CL>=.975,1,'First') length(CL)]));
    AC.IS_SCNA(i,1)=(AC.CCF_Low(i,1)>0)&(abs(AC.tau(i)-NCN)>4*AC.sigma_tau(i));
    AC.CCF_prob_nextbest(i,:)=Lnext/sum(Lnext);
    if ~isnan(NA2+NB2)  
       AC.CCF_hat_nextbest(i,1)=CCF(find(squeeze(chi(NA2+1,NB2+1,:))==min(squeeze(chi(NA2+1,NB2+1,:))),1,'first'));
    end
    AC.purity(i,1)=purity;
    AC.ploidy(i,1)=ploidy;
    AC.muA_fit(i,1)=muA1;
    AC.muB_fit(i,1)=muB1;
    AC.tau_fit(i,1)=tau1;
    AC.muA_nextbest(i,1)=muA2;
    AC.muB_nextbest(i,1)=muB2;
    AC.tau_nextbest(i,1)=tau2;
    AC.NCN(i,1)=NCN;
    
    % somatic something happened? 
    af_balance_signal=abs(AC.mu_major(i)-1)/AC.sigma_major(i);
    if isnan(af_balance_signal)||(NCN<2)
        af_balance_signal=0;
    end    
    scna_signal=abs(AC.tau(i)-NCN)/AC.sigma_tau(i);
    AC.IS_SCNA(i,1)=(AC.CCF_Low(i,1)>0)& ((scna_signal>PAR.SIGMA_THRESHOLD)|(af_balance_signal>PAR.SIGMA_THRESHOLD));
    i;
end

X=AC;

end

%% 


function test
%% PCAWG
AC1=load_table('/local/cga-fh/cga/PCAWG_train2_pipeline_B/Pair/RK175_C01_RK175_B01/jobs/capture/AllelicCapseg/latest/results/RK175_C01_RK175_B01.tsv')
purity1=0.636; ploidy1=1.966; gender1='M'
PAR=[];
PAR.CCF_BINS=0:0.01:1;
PAR.SMOOTH_SPAN=7;
PAR.HIGH_CCF_CHISQ=1;
PAR.CN_CHISQ=0.5;
PAR.MAXCN=15;
PAR.SIGMA_MIN=0.05;
PAR.SIGMA_THRESHOLD=2;
PAR.X_MALE_TAU_SCALE = 2
X1=AllelicCapSeg_PP_CCF_fit(AC1,purity1,ploidy1,gender1,PAR)

%% DLBCL 

clear
AC=load_tsv('/Users/Stewart/Projects/Cancer/DLBCL/ABSOLUTE/DLBCL_290_22Mar2015.28Mar2015.allelic_capseg.tsv')
P=load_tsv('/Users/Stewart/Projects/Cancer/DLBCL/FH/DLBCL_290_22Mar2015.28Mar2015.pairs.tsv')
T=load_tsv('/Users/Stewart/Projects/Cancer/DLBCL/FH/DLBCL_290_22Mar2015.28Mar2015.segtab.tsv')
G=load_tsv('/Users/stewart/Projects/Cancer/DLBCL/FH/gender_379.15Mar2015.txt')
P.sample=cellfun(@(x) x(1),regexp(P.pair_id,'-TP-','split'))
P.sample=regexprep(P.sample,'-nullpair$','')
[i m]=ismember(AC.pair_id,P.pair_id)
AC.sample(i,1)=P.sample(m(i))
a=tab(AC.sample)
s=a.x{a.n==max(a.n)};
s=a.x{a.n==min(a.n)};
AC1=trimStruct(AC,ismember(AC.sample,s))
purity1=P.purity(ismember(P.sample,s))
ploidy1=P.ploidy(ismember(P.sample,s))
gender1=G.sex{ismember(G.sample,s)}

X1=AllelicCapSeg_PP_CCF_fit(AC1,purity1,ploidy1,gender1)

tab(X1.NA*100+X1.NB)
k=find(X1.NA==2 & X1.NB==3)
k=find(X1.NA==2 & X1.NB==2)
[purity1, ploidy1, ismember(gender1,'F')]
myhist(X1.tau,0:0.1:5,'log')
hold on; myhist(X1.tau(k),0:0.1:5,'log'); hold off

plot(X1.tau,X1.tau_fit,'+',X1.tau(k),X1.tau_fit(k),'rx')
plot(X1.mu_major,X1.muB_fit,'+',X1.mu_major(k),X1.muB_fit(k),'rx')

myhist(X1.CCF_hat,0:0.01:1,'log')
hold on; myhist(X1.CCF_hat(k),0:0.01:1,'log'); hold off

tab(X1.IS_SCNA)
tab(X1.IS_SCNA(k))

myhist(X1.mu_minor,0:0.1:5,'log')
hold on; myhist(X1.mu_minor(k),0:0.1:5,'log'); hold off
myhist(X1.mu_major,0:0.1:5,'log')
hold on; myhist(X1.mu_major(k),0:0.1:5,'log'); hold off
myhist(X1.mu_major,0:0.1:5,'log')
hold on; myhist(X1.mu_major(k),0:0.1:5,'log'); hold off

s

round([ 100*purity1 100*ploidy1 length(X1.IS_SCNA) 100*mean(X1.IS_SCNA) 100*sum(X1.IS_SCNA.*X1.length)/sum(X1.length) sum((X1.NA==2)&(X1.NB==2))])
xL=xhg19(1:25,0*(1:25))'/1e6;

z1=X1;
z1.x1=xhg19(z1.Chromosome,z1.Start_bp)/1e6;
z1.x2=xhg19(z1.Chromosome,z1.End_bp)/1e6;
z2=trimStruct(z1,z1.IS_SCNA);
clf
subplot(3,1,1)
semilogy([z1.x1 z1.x2]',[z1.tau z1.tau]','b-',[z2.x1 z2.x2]',[z2.tau z2.tau]','r-','linewidth',2)
line([xL xL]',[0.05+0*xL 20+0*xL]','color',0.25*[1 1 1],'linestyle','--')
set(gca,'ylim',[0.05 20],'ytick',[0.1 1 10],'xlim',[0 3000])
ylabel('CN ratio')
title(sprintf('%s  %s  purity=%.2f ploidy=%.2f',z1.sample{1},gender1,purity1,ploidy1) );
subplot(3,1,2)
semilogy([z1.x1 z1.x2]',[z1.mu_major z1.mu_major]','b-',[z2.x1 z2.x2]',[z2.mu_major z2.mu_major]','r-',...
    [z1.x1 z1.x2]',[z1.mu_minor z1.mu_minor]','b-',[z2.x1 z2.x2]',[z2.mu_minor z2.mu_minor]','r-','linewidth',2)
line([xL xL]',[0.05+0*xL 20+0*xL]','color',0.25*[1 1 1],'linestyle','--')
set(gca,'ylim',[0.05 20],'ytick',[0.1 1 10],'xlim',[0 3000])
ylabel('Allele fraction ratios')
subplot(3,1,3)
yA1=0.05+[z1.NA.*z1.CCF_hat+1*(1-z1.CCF_hat) z1.NA.*z1.CCF_hat+1*(1-z1.CCF_hat)]';
yA2=0.05+[z2.NA.*z2.CCF_hat+1*(1-z2.CCF_hat) z2.NA.*z2.CCF_hat+1*(1-z2.CCF_hat)]';
yB1=0.05+[z1.NB.*z1.CCF_hat+1*(1-z1.CCF_hat) z1.NB.*z1.CCF_hat+1*(1-z1.CCF_hat)]';
yB2=0.05+[z2.NB.*z2.CCF_hat+1*(1-z2.CCF_hat) z2.NB.*z2.CCF_hat+1*(1-z2.CCF_hat)]';
semilogy([z1.x1 z1.x2]',yA1,'b-',[z2.x1 z2.x2]',yA2,'r-',...
    [z1.x1 z1.x2]',yB1,'b-',[z2.x1 z2.x2]',yB2,'r-','linewidth',2)
line([xL xL]',[0.05+0*xL 20+0*xL]','color',0.25*[1 1 1],'linestyle','--')
set(gca,'ylim',[0.05 20],'ytick',[0.1 1 10],'xlim',[0 3000])
ylabel('fit Allele CN')


end
