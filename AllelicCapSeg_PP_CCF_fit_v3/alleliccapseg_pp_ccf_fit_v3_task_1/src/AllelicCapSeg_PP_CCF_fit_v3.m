function [X]=AllelicCapSeg_PP_CCF_fit_v3(AC,purity,ploidy,gender,PAR)
%[X]=AllelicCapSeg_PP_CCF_fit_v3 (AC,purity,ploidy,PAR)
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
%       NA: SCNA Allele A copy number (lower or equal aCN)
%       NB: SCNA Allele B copy number (higher or equal aCN)
%       chisq: fit chisquared
%       CCF_NA: CCF estimate NA
%       CCF_NB: CCF estimate NA
%       CCF_NA_Low: CCF_NA 95% low CI
%       CCF_NA_High: CCF_NA 95% high CI
%       CCF_NB_Low: CCF_NB 95% low CI
%       CCF_NB_High: CCF_NB 95% high CI
%       IS_SCNA: tumor SCNA flag
%       purity: 
%       ploidy: 
%       gender: MALE or FEMALE
%       CNN: normal CN state (2, except 1 for male X or Y, 0 for female Y)
%       muA_fit: fit value of mu_minor
%       muB_fit: fit value of mu_major
%       tau_fit: fit value of tau
%
%  Chip 14 Mar 20
%
if ~isstruct(AC)
    AC=load_tsv(AC)
end
x1=xhg19(AC.Chromosome,AC.Start_bp);
[~,k]=sort(x1);
AC=trimStruct(AC,k);

if (nargin<5)
    PAR=[];
    
    PAR.ACN_BINWIDTH=0.0025;
    PAR.ACN_MAX=7;
    PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD=0.9
    PAR.SIGMA_MIN=0.025;
    PAR.X_MALE_TAU_SCALE = 1
    PAR.min_ACNRES=0.05
end
if (nargin<4)
    msgID = 'AllelicCapSeg_PP_CCF_fit:Input';
    msg = 'Not enough inputs: AC, purity, ploidy, gender';
    throw(MException(msgID,msg));
end
if ~isfield(PAR,'ACN_BINWIDTH')
    PAR.ACN_BINWIDTH=0.0025;
end
if ~isfield(PAR,'ACN_MAX')
    PAR.ACN_MAX=7;
end
if ~isfield(PAR,'JOINT_OVERLAP_CLUSTERING_THRESHOLD')
    PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD=0.9
end

if ~isfield(PAR,'SIGMA_MIN')
    PAR.SIGMA_MIN=0.025; 
end
if ~isfield(PAR,'SIGMA_THRESHOLD')
     PAR.SIGMA_THRESHOLD=2;
end
if ~isfield(PAR,'min_ACNRES')
    PAR.min_ACNRES=0.05
end
    
ACRBINS=0:PAR.ACN_BINWIDTH:PAR.ACN_MAX;
ACRCONV=-0.5:PAR.ACN_BINWIDTH:0.5;


% set mu for nan segments with sigma_tau<diff(tau) 
% assume missing mu's are simply noise+missing hets w/o aCN change
k=find(isnan(AC.mu_minor))
d1=1e6+0*k
k1=k(k>1)
i1=AC.Chromosome(k1)==AC.Chromosome(k1-1)
k1=k1(find(i1))
d1(find(i1))=abs(AC.tau(k1)-AC.tau(k1-1))
d2=1e6+0*k
k2=k(k<length(AC.mu_minor))
i2=AC.Chromosome(k2)==AC.Chromosome(k2+1)
k2=k2(find(i2))
d2(find(i2))=abs(AC.tau(k2)-AC.tau(k2+1))
[d iside]=min([d1 d2],[],2)
iside(iside==1)=-1
iside(iside==2)=1
ks=find(d<AC.sigma_tau(k))
AC.mu_minor(k(ks))
AC.mu_minor(k(ks))=AC.mu_minor(k(ks)+iside(ks))
AC.mu_major(k(ks))=AC.mu_major(k(ks)+iside(ks))
AC.sigma_minor(k(ks))=AC.sigma_tau(k(ks))
AC.sigma_major(k(ks))=AC.sigma_tau(k(ks))

% set nan sigmas to something other than nan
k=find(isnan(AC.sigma_major))
AC.sigma_minor(k)=nanmedian(AC.sigma_minor)
AC.sigma_major(k)=nanmedian(AC.sigma_major)

% add minimal sigma to each segment
AC.sigma_major=sqrt(AC.sigma_major.^2+PAR.SIGMA_MIN.^2);
AC.sigma_minor=sqrt(AC.sigma_minor.^2+PAR.SIGMA_MIN.^2);
AC.sigma_tau=sqrt(AC.sigma_tau.^2+PAR.SIGMA_MIN.^2);

% intitialize AC fields
AC.NA=NaN*AC.length;
AC.NB=NaN*AC.length;
% AC.chisq=NaN*AC.length;
AC.NA_CCF=NaN*AC.length;
% AC.CCF_NA_Low=NaN*AC.length;
% AC.CCF_NA_High=NaN*AC.length;
AC.NB_CCF=NaN*AC.length;
% AC.CCF_NB_Low=NaN*AC.length;
% AC.CCF_NB_High=NaN*AC.length;
% AC.muA_fit=NaN*AC.length;
% AC.muB_fit=NaN*AC.length;
% AC.tau_fit=NaN*AC.length;
AC.purity=purity+0*AC.length;
AC.ploidy=ploidy+0*AC.length;
AC.gender=repmat(gender,length(AC.length),1);
AC.CCN=2+0*AC.length;
AC.IS_SCNA=AC.length<-100;

XCN=2; % X normal CN 
if ismember(upper(gender),{'M','MALE'})
    XCN=1;    
    k=find(AC.Chromosome>=23)
    AC.tau(k)=AC.tau(k)/PAR.X_MALE_TAU_SCALE;
    AC.mu_minor(k)=0;
    AC.mu_major(k)=AC.tau(k)
    AC.CCN(k)=1
    
end

% ACR bin for each segment 
[nA,bA]=histc(AC.mu_minor,ACRBINS);
[nA,bB]=histc(AC.mu_major,ACRBINS);
[nT,bT]=histc(AC.tau,ACRBINS);
N=length(AC.tau);
x=zeros(N,length(ACRBINS));w=x; 
xA=x; wA=x; xB=x; wB=x; xT=x; wT=x;
% loop over segments 
for i1=1:N
    if (bA(i1)>0)
        xA(i1,bA(i1))=AC.length(i1)/1e6;
        xB(i1,bB(i1))=AC.length(i1)/1e6;
        x=xA+xB;
        sig=normpdf(ACRCONV,0,AC.sigma_minor(i1)); sig=sig/sum(sig);
        wA(i1,:) = conv(xA(i1,:),sig,'same');
        wB(i1,:) = conv(xB(i1,:),sig,'same');
        w=wA+wB;
    end
%     xT(i1,bT(i1))=AC.length(i1)/1e6;
%     sig=normpdf(ACRCONV,0,AC.sigma_tau(i1)); sig=sig/sum(sig);
%     wT(i1,:) = conv(xT(i1,:),sig,'same');          
end

% marginal segment CR sum
xs=sum(x,1);
% marginal segment CR sum smoothened by sigma
ws=sum(w,1);
% mean allelic copy ratio
avecr=sum(ACRBINS.*ws)/sum(ws);
% total CR
wt=sum(wT,1);
% mean copy ratio
% aveTCR=sum(ACRBINS.*wt)/sum(wt);

% distance between integer copy levels
dcr=avecr/(ploidy/2+(1-purity)/purity);
% aCR for copy number 0
cr0=dcr*(1-purity)/purity;
% map ACRBINS to ACNBINS
ACNBINS=(ACRBINS-cr0)/dcr;
% allelic copy ratios of integer somatic copy numbers 
ACR=(0:PAR.ACN_MAX)*dcr+cr0

AC.cn_minor=(AC.mu_minor-cr0)/dcr;
AC.cn_major=(AC.mu_major-cr0)/dcr;
AC.cn_sigma_minor=AC.sigma_minor/dcr;
AC.cn_sigma_major=AC.sigma_major/dcr;

% can we adjust the conversion (dcr, cr0) such that a minimum number of 
% bases are just above integer CN levels? 

% plot(x1,AC.cn_minor,'o',x1,AC.cn_major,'o')
% errorbar(x1,AC.cn_minor,AC.cn_sigma_minor,'.'); hold on; errorbar(x1,AC.cn_major,AC.cn_sigma_major,'.');hold off


[ga, gb]=CNcluster(AC,wA,wB,PAR);
%
GA=unique(ga);
for n=1:length(GA)
    k=find(ga==GA(n));
    AC.cn_clusterA(k,1)=GA(n);
    w=AC.length(k); w=w/sum(w);
    cn=sum(AC.cn_minor(k).*w);
    AC.cn_minor_cluster(k,1)=cn;
    AC.cn_minor_sigma_cluster(k,1) = sqrt(var(AC.cn_minor(k),w)+sum(AC.cn_sigma_minor(k).*w).^2);
end
GB=unique(gb);
for n=1:length(GB)
    k=find(gb==GB(n));
    AC.cn_clusterB(k,1)=GB(n);
    w=AC.length(k); w=w/sum(w);
    cn=sum(AC.cn_major(k).*w);
    AC.cn_major_cluster(k,1)=cn;
    AC.cn_major_sigma_cluster(k,1) = sqrt(var(AC.cn_major(k),w)+sum(AC.cn_sigma_major(k).*w).^2);    
end



%   for g=GA'
%       k=find(AC.cn_GA==g);
%       [g length(k)]
%       [AC.cn_minor_cluster(k) AC.cn_minor(k) AC.mu_minor(k) AC.length(k)/1e6]
%   end
%   
%   for g=GB'
%       k=find(AC.cn_GB==g);
%       [g length(k)]
%       [AC.cn_major_cluster(k) AC.cn_major(k) AC.mu_major(k) AC.length(k)/1e6]
%   end
%   
%   WIN=0.1
%   options = optimset('TolX',1e-6)% ,'Display','iter')%'optimplotfval','@optimplotfval')
%   cA=AC.cn_minor_cluster; cB=AC.cn_major_cluster;  w=AC.length/1e6;
%   cA=AC.cn_minor; cB=AC.cn_major;  w=AC.length/1e6;
%   fun = @(x)cost0(x,cA,cB,w,WIN);
%   p0 = [0.99];
%   [p,fval,exitflag,output] = fminsearch(fun,p0,options)
%   histw([AC.length; AC.length],[p*AC.cn_minor; p*AC.cn_major],-0.1:0.01:4)
%   histw([AC.length; AC.length],[AC.cn_minor; AC.cn_major],-0.1:0.01:4)

% minor shift in ploidy such that somatic copy numbers do not fluctuate over integers  
W=0.1;
cA=AC.cn_minor_cluster; cB=AC.cn_major_cluster;  len=AC.length;
fun = @(x)cost1(x,cA,cB,len,W);
p0 = [0.99];
A = fminsearch(fun,p0);

AC.cn_minor_cluster=A*AC.cn_minor_cluster;
AC.cn_major_cluster=A*AC.cn_major_cluster;
AC.cn_minor=A*AC.cn_minor;
AC.cn_major=A*AC.cn_major;
AC.cn_ploidy=AC.ploidy/A;
DCR=A*dcr;

% fix cn_minor_cluster that fall just below 0
k=find( (AC.cn_minor_cluster<0)&(AC.cn_minor_cluster>(-2*AC.cn_minor_sigma_cluster)));
AC.cn_minor_cluster(k)=0;
  
N=length(AC.purity);
for i=1:N
    
    aCNA=AC.cn_minor_cluster(i);
    aCNAW=sqrt(PAR.min_ACNRES^2+((AC.sigma_minor(i))/DCR)^2);

    b=sort([-10 (0:10)-aCNAW (0:10)+aCNAW]);
    bCN=[NaN sort([0 0 0 2:10 2:10])];
    
    [n1 ibin]=histc(aCNA,b);
    switch ibin
        case 2  
            AC.NA(i)=0;
            AC.NA_CCF(i)=1;
            AC.IS_SCNA(i)=1>0;
        case 3 
            AC.NA(i)=0;
            AC.NA_CCF(i)=1-aCNA;
            AC.IS_SCNA(i)=1>0;
        case 4 
            AC.NA(i)=0;
            AC.NA_CCF(i)=0;
        case 5 
            AC.NA(i)=2;
            AC.NA_CCF(i)=aCNA-1;
            AC.IS_SCNA(i)=1>0;
        case {6,8,10,12,14,16,18,20}
            AC.NA(i)=round(aCNA);
            AC.NA_CCF(i)=1;
            AC.IS_SCNA(i)=1>0;
        case {7,9,11,13,15,17,19,21}
            AC.NA(i)=ceil(aCNA);
            AC.NA_CCF(i)=(aCNA-1)/(AC.NA(i)-1);
            AC.IS_SCNA(i)=1>0;
            
    end
    
    aCNB=AC.cn_major_cluster(i);
    aCNAW=sqrt(PAR.min_ACNRES^2+((AC.sigma_minor(i))/DCR)^2);

    b=sort([-10 (0:10)-aCNAW (0:10)+aCNAW]);
    bCN=[NaN sort([0 0 0 2:10 2:10])];
    
    [n1 ibin]=histc(aCNB,b);
    switch ibin
        case 2  
            AC.NB(i)=0;
            AC.NB_CCF(i)=1;
            AC.IS_SCNA(i)=1>0;
        case 3 
            AC.NB(i)=0;
            AC.NB_CCF(i)=1-aCNB;
            AC.IS_SCNA(i)=1>0;
        case 4 
            AC.NB(i)=0;
            AC.NB_CCF(i)=0;
        case 5 
            AC.NB(i)=2;
            AC.NB_CCF(i)=aCNB-1;
            AC.IS_SCNA(i)=1>0;
        case {6,8,10,12,14,16,18,20}
            AC.NB(i)=round(aCNB);
            AC.NB_CCF(i)=1;
            AC.IS_SCNA(i)=1>0;
        case {7,9,11,13,15,17,19,21}
            AC.NB(i)=ceil(aCNB);
            AC.NB_CCF(i)=(aCNB-1)/(AC.NB(i)-1);
            AC.IS_SCNA(i)=1>0;
            
    end
    
end

AC.DeltaCR=dcr+0*AC.purity;
AC.averageCR=avecr+0*AC.purity;
AC.zeroCR=cr0+0*AC.purity;

X=AC;

if(0)
    A
    subplot(3,1,1)
    histw([AC.length; AC.length],[AC.cn_minor; AC.cn_major],-0.1:0.01:4)
    subplot(3,1,2)
    histw([AC.length; AC.length],[A*AC.cn_minor_cluster; A*AC.cn_major_cluster],-0.1:0.01:4)
    subplot(3,1,3)
    histw([AC.length; AC.length],[A*AC.cn_minor; A*AC.cn_major],-0.1:0.01:4)
end

end


function [c]=cost0(p,cA,cB,len,W)
L=sum(len);
c=0;
a=p(1);
for CN=0:7
    % small window above integer CN
    kA=find(((cA*a)>(CN-W))&((cA*a)<(CN+W)));
    kB=find(((cB*a)>(CN-W))&((cB*a)<(CN+W)));
    k=unique([kA; kB]);
    wA=1*ismember(k,kA);
    wB=1*ismember(k,kB);
    c = c - sum( (len(k)/L).*(wA.*(cA(k)-CN)).^2+(wB.*(cB(k)-CN).^2));
end
c = c + 0.0001*(a-1).^2;
[p c]
end


function [c]=cost1(p,cA,cB,len,W)
L=sum(len);
c=0;
a=p(1);

for CN=2:5 % 0, 1 ok
    % small window above integer CN
    kA=find(((cA*a)>CN)&((cA*a)<(CN+W)));
    kB=find(((cB*a)>CN)&((cB*a)<(CN+W)));
    k=unique([kA; kB]);
    c = c + (sum(len(k))/L).^2; 
end
c = c + ((0.001+100*(a>1))*(a-1)).^2 ;
end

function [c]=cost2(p,cA,cB,len,W)
L=sum(len);
c=0;
a=p(1);
b=p(2);

kA=find(((cA*a+b)<0)&((cA*a+b)>-W));
kB=find(((cB*a+b)<0)&((cB*a+b)>-W));
k=unique([kA; kB]);
c = c + (sum(len(k))/L).^2;
for CN=2:5 % 0, 1 ok
    % small window above integer CN
    kA=find(((cA*a+b)>CN)&((cA*a+b)<(CN+W)));
    kB=find(((cB*a+b)>CN)&((cB*a+b)<(CN+W)));
    k=unique([kA; kB]);
    c = c + (sum(len(k))/L).^2;
end
c = c + (a-1).^2 + b.^2;
end


function [GA, GB]=CNcluster(AC,wA,wB,PAR)
if (nargin<5)
    PAR=[];
    PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD=0.9;
end
if ~isfield(PAR,'JOINT_OVERLAP_CLUSTERING_THRESHOLD')
    PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD=0.9;
end
N=length(AC.cn_minor)
[~,kA]=sort(AC.cn_minor);
[~,kB]=sort(AC.cn_major);
%  plot(x1,AC.cn_minor(kA),'o')
classA=1:N;
nxtA=classA; 
classB=classA; nxtB=nxtA; % self-linked list
for i=1:(N-1)
    k1=kA(i);
    k2=kA(i+1);
    w1=wA(k1,:); w1=w1/sum(w1);
    w2=wA(k2,:); w2=w2/sum(w2);
    w12=(w1.*w2);
    x1=sum(w12)/sum(w1.^2);
    x2=sum(w12)/sum(w2.^2);
    % plot(w1); hold on; plot(w2); plot(w12); hold off; xlim([100 400]); [x1 x2]
    if (max(x1,x2)>PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD)
        [classA,nxtA]=connect(k1,k2,classA,nxtA);
    end 
    k1=kB(i);
    k2=kB(i+1);
    w1=wB(k1,:); w1=w1/sum(w1);
    w2=wB(k2,:); w2=w2/sum(w2);
    w12=(w1.*w2);
    x1=sum(w12)/sum(w1.^2);
    x2=sum(w12)/sum(w2.^2);
    % plot(w1); hold on; plot(w2); plot(w12); hold off; xlim([100 400]); [x1 x2]
    if (max(x1,x2)>PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD)
        [classB,nxtB]=connect(k1,k2,classB,nxtB);
    end       
    
end
UG=unique(classA);
GA=zeros(N,1);
for n=1:length(UG)
    k=find(classA==UG(n));
    GA(k)=n+0*k;
end
UG=unique(classB);
GB=zeros(N,1);
for n=1:length(UG)
    k=find(classB==UG(n));
    GB(k)=n+0*k;    
end
end

function [class,nxt]=connect(i,j,class,nxt)
% Knuth Equivalence Class algorithm
if (class(i)==class(j)), return;end
j1=j;
class(j1)=class(i);
while (nxt(j1)~=j)
    j1=nxt(j1);
    class(j1)=class(i);
end
i1=nxt(i);
nxt(i)=j;
nxt(j1)=i1;
end

function test
%%
clear
clf
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/CN/test')

ACF='REBC-ACCF-TTP1-A-1-1-D-A49W-36.modelFinal.seg.final.seg.pruned.seg.merged.seg.acs.seg' ; purity=0.28	; ploidy=2.49; gender='male'
ACF='REBC-ACAB-TTP1-A-1-1-D-A49V-36.modelFinal.seg.final.seg.pruned.seg.merged.seg.acs.seg' ; purity=0.58; ploidy=5.45; gender='male'
ACF='REBC-AC9P-TTP1-A-1-1-D-A49U-36.modelFinal.seg.final.seg.pruned.seg.merged.seg.acs.seg' ; purity=0.792; ploidy=2; gender='male'
ACF='REBC-AC9O-TTP1-A-1-1-D-A49U-36.modelFinal.seg.final.seg.pruned.seg.merged.seg.acs.seg' ; purity=0.688; ploidy=2; gender='male'
ACF='REBC-AF72-TTP1-A-1-1-D-A649-36.modelFinal.seg.final.seg.pruned.seg.merged.seg.acs.seg';purity=0.844	; ploidy=1.75; gender='male'; purity=0.42	; ploidy=1.81;  %; ploidy=1.75;  



AC=load_tsv(ACF)
PAR=[]
PAR.JOINT_OVERLAP_CLUSTERING_THRESHOLD=0.01
X1=AllelicCapSeg_PP_CCF_fit_v3(AC,purity,ploidy,gender)

end

function test2
%%
clear
clf
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/CN')
A=load_tsv('REBC_primary_pair_393.cnv_postprocessing_tumor_acs.tsv')
P=load_tsv('/Volumes/GoogleDrive/My Drive/Cancer/REBC/FC/REBC.pairs.15Mar2020.tsv')

PT=load_tsv('/Volumes/GoogleDrive/My Drive/Cancer/REBC/FC/REBC.participants.tsv')
[i m]=ismember(P.participant,PT.entity_participant_id);
P.Gender(i,1)=PT.Gender(m(i));
pid1=unique(sort(A.sample));

[i m]=ismember(P.entity_pair_id,pid1);
P=trimStruct(P,i)
N=length(pid1)
AC=[]
for i=1:N

    A1=trimStruct(A,ismember(A.sample,pid1(i)));
    P1=trimStruct(P,ismember(P.entity_pair_id,pid1(i)))    ;
    X1=AllelicCapSeg_PP_CCF_fit_v3(A1,P1.purity(1),P1.absolute_ploidy(1),P1.Gender(1));
    if isempty(AC)
       AC=X1;
    else
       AC=mergeStruct(AC,X1);
    end
    [i X1.N]
end    
printStruct(AC,-1,'CCF/REBC_primary_pair_393.cnv_postprocessing_tumor_acs.ccf.tsv')


end


