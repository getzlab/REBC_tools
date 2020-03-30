function [ X1, XA ] = MAF_AC_PP_CCF_fit_v3(X,T,purity,ploidy,gender,PAR)
%  [ X1 , XA ] = MAF_AC_PP_CCF_fit_v3(X, T ,purity,ploidy,gender,PAR)
%
%  MAF_AC_PP_CCF_fit_v3 adds estimated mutation CCF and local tumor SCNA allele copy number
%  [NA,NB] and  CCF's to AC based on puity and ploidy (from ABSOLUTE or other method).
%  This function makes estimates for all segments, even X and Y when gender
%  info is available from PAR.
%  MAF_AC_PP_CCF_fit_v3 operates on one sample at a time (one purity, ploidy, gender)
%
%  Note: X1 output roughly correpsonds to the ABSOLUTE ABS_MAF file, with a
%  few caveats:
%      1) mutation multiplicity is 1 or NA or NB (no intermediate m) with a
%      prior that as the SCNA CCF goes up the likelihood that the SCNA preceeded
%      the mutation increases.
%      2) When allele fraction betapdf's exceed CCF=1, the overflow is assigned to CCF=1
%
%  XA includes ABSOLUTE maf file fields - with caveats (above) and some
%  blank fields (eg. Pr_GL_som_HZ_alt - is there a document about ABSOLUTE maf fields?) 
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
if ischar(X)
    X=load_table(X);
end
if ischar(T)
    T=load_table(T);
end
if nargin<6
    PAR=[];
    PAR.DEFAULT=1;
end
if ~isfield(PAR,'DELTA')
    PAR.DELTA=0.01;
end
if ~isfield(PAR,'OUTPUT_MAF')
    PAR.OUTPUT_MAF='';
end
if ~isfield(PAR,'CONVERT_PCCF_TO_FIELDS')
    PAR.CONVERT_PCCF_TO_FIELDS=0;
end
if ~isfield(PAR,'SMOOTH')
    PAR.SMOOTH=50001;
end
if ~isfield(PAR,'KEEP_PCCF_VECTOR')
    PAR.KEEP_PCCF_VECTOR=1;
end
if ~isfield(PAR,'ASSUME_NORMAL_X_WHEN_MISSING')
    PAR.ASSUME_NORMAL_X_WHEN_MISSING=1;
end
if ~isfield(PAR,'CCF_PRIOR')
    PAR.CCF_PRIOR_SLOPE=0.01;
end
if ~isfield(PAR,'ADD_ABSOLUTE_FIELDS')
    PAR.ADD_ABSOLUTE_FIELDS=1;
end

XA=[];

AF_BINS=0:PAR.DELTA:1;
CCF_BINS=0:PAR.DELTA:4; % allow for fluctuations beyond CCF=1 (trim later)

XCNA=1; % X normal CN  A allele
XCNB=1; % X normal CN  B allele
if ismember(upper(gender),{'M','MALE'})
    XCNB=0; % no B allele
end


X.SCNA_NA=NaN*X.Start_position;
X.SCNA_NB=NaN*X.Start_position;
X.SCNA_q_hat=NaN*X.Start_position;
X.SCNA_NA_CCF=NaN*X.Start_position;
X.SCNA_NA_CCF_Low=NaN*X.Start_position;
X.SCNA_NA_CCF_High=NaN*X.Start_position;
X.SCNA_NB_CCF=NaN*X.Start_position;
X.SCNA_NB_CCF_Low=NaN*X.Start_position;
X.SCNA_NB_CCF_High=NaN*X.Start_position;
X.SCNA_tau=NaN*X.Start_position;
X.SCNA_NCNA=NaN*X.Start_position;
X.SCNA_NCNB=NaN*X.Start_position;
X.IS_SCNA=NaN*X.Start_position;

X.purity=purity*ones(size(X.Start_position));
X.ploidy=ploidy*ones(size(X.Start_position));

X.dna_fraction_in_tumor=NaN*X.Start_position;
X.ccf_hat=NaN*X.Start_position;
X.ccf_CI_low=NaN*X.Start_position;
X.ccf_CI_high=NaN*X.Start_position;

X.ccf_CI95_low=NaN*X.Start_position;
X.ccf_CI95_high=NaN*X.Start_position;
X.ccf_mode=NaN*X.Start_position;
X.ccf_mean=NaN*X.Start_position;
X.ccf_median=NaN*X.Start_position;
X.pccf=NaN(length(X.Start_position),length(AF_BINS));
X.clonal = NaN*X.Start_position;
X.multiplicity = NaN*X.Start_position;
X.ccf_ALT_model = NaN*X.Start_position;
X.ccf_ALT_model_AF = NaN*X.Start_position;
X.ccf_ALT_model_OVER = NaN*X.Start_position;

% % CCF_PRIOR (prefer higher CCFs from max likelihood to break ties)
CCF_PRIOR=1+0*CCF_BINS; CCF_PRIOR(CCF_BINS>1)=0;  CCF_PRIOR=CCF_PRIOR./sum(CCF_PRIOR);
CCF_PRIOR=(CCF_BINS-0.5)*PAR.CCF_PRIOR_SLOPE*CCF_PRIOR(1)+CCF_PRIOR;
CCF_PRIOR=CCF_PRIOR/mean(CCF_PRIOR(CCF_BINS<=1));CCF_PRIOR(CCF_BINS>1)=0;

% alt allele fraction
if ~isfield(X,'i_tumor_f')
    X.i_tumor_f=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
end


if isnumeric(X.Chromosome)
    X.Chromosome=num2chrom(X.Chromosome);
end

N=length(X.Start_position);

T.NCNA=1+0*T.Start_bp;
T.NCNB=1+0*T.Start_bp;

for i=1:N
    
    NCNA=XCNA;
    
    C=X.Chromosome(i);
    if ~isnumeric(X.Chromosome)
        C=chrom2num(C);
    end
    
    k=find(ismember(T.Chromosome,C)&(T.Start_bp<=X.Start_position(i))&(T.End_bp>=X.Start_position(i)));
    if numel(k)==0
        % no SCNA segment
        if (C>22)&&PAR.ASSUME_NORMAL_X_WHEN_MISSING
            % X or Y
            NCNB=XCNB;
            if (C==24)&&ismember(upper(gender),{'F','FEMALE'})
                %NCNA=0;
                continue % no female Y
            end
            k=length(T.NA)+1;
            T1=trimStruct(T,1);
            T1.Chromosome=23;
            [~, x2]=xhg19(xhg19(24,1)-2);
            T1.Start_bp=1;
            T1.End_bp=x2;
            T1.tau=NaN;
            T1.mu_minor=NaN;
            T1.mu_major=NaN;
            T1.cn_minor=NaN;
            T1.cn_major=NaN;
            T1.cn_clusterA=NaN;
            T1.cn_clusterB=NaN;
            T1.cn_minor_sigma_cluster=0;
            T1.cn_major_sigma_cluster=0;
            T1.NCNA=NCNA;
            T1.NCNB=NCNB;
            T1.NA=0;
            T1.NB=0;
            T1.NA_CCF=0;
            T1.NB_CCF=0;
            T1.IS_SCNA=NaN;
            T=mergeStruct(T,T1);
            
        else
            continue;
        end
    end
    % if two segments overlap (problem) pick 1
    if numel(k)>1, k=k(1); end
    % assign SCNA state at mutation site
    X.SCNA_NA(i)=T.NA(k);
    X.SCNA_NB(i)=T.NB(k);
    X.SCNA_NA_CCF(i)=T.NA_CCF(k);
    X.SCNA_NA_CCF_Low(i)=max(0,T.NA_CCF(k)-T.cn_minor_sigma_cluster(k));
    X.SCNA_NA_CCF_High(i)=min(1,T.NA_CCF(k)+T.cn_minor_sigma_cluster(k));
    X.SCNA_NB_CCF(i)=T.NB_CCF(k);
    X.SCNA_NB_CCF_Low(i)=max(0,T.NB_CCF(k)-T.cn_major_sigma_cluster(k));
    X.SCNA_NB_CCF_High(i)=min(1,T.NB_CCF(k)+T.cn_major_sigma_cluster(k));
    X.SCNA_tau(i)=T.tau(k);
    X.IS_SCNA(i)=T.IS_SCNA(k);
    X.SCNA_NCNA(i)=T.NCNA(k);
    X.SCNA_NCNB(i)=T.NCNB(k);
    X.SCNA_q_hat(i)=X.SCNA_NA(i)+X.SCNA_NB(i);
    
    % DNA per cell
    tumorDNA=X.SCNA_NA(i).*X.SCNA_NA_CCF(i)+X.SCNA_NB(i).*X.SCNA_NB_CCF(i).*X.SCNA_NCNB(i)+X.SCNA_NCNA(i).*(1-X.SCNA_NA_CCF(i))+X.SCNA_NCNB(i).*(1-X.SCNA_NB_CCF(i));
    if tumorDNA<=0
        %no tumor DNA -> no mutation
        fprintf('tumorDNA=0 chr:%s pos:%d NA=%d NA_CCF=%.3f NB=%d NB_CCF-%.3f ALT=%d REF=%d',X.Chromosome{i},X.Start_position(i),X.SCNA_NB(i),X.SCNA_NA_CCF(i),X.SCNA_NB(i),X.SCNA_NB_CCF(i),X.t_alt_count(i),X.t_ref_count(i))
        X.ccf_Low(i,1)=0;
        X.ccf_High(i,1)=0;
        X.ccf_CI95_low(i,1)=0;
        X.ccf_CI95_high(i,1)=0;
        X.ccf_mode(i,1)=0;
        X.ccf_mean(i,1)=0;
        X.ccf_median(i,1)=0;
        X.multiplicity(i,1)=0;
        X.ccf_hat(i,1)=0;
        X.clonal(i,1) = 0;
        continue
    end
    normalDNA=X.SCNA_NCNA(i)+X.SCNA_NCNB(i);
    totalDNA = tumorDNA*purity + normalDNA*(1-purity);
    
    X.dna_fraction_in_tumor(i)=(purity*tumorDNA)./(purity*tumorDNA+(1-purity)*normalDNA);
    
    % skip events that ABSOLUTE skipped
    if isnan(X.dna_fraction_in_tumor(i)), continue; end
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
    ccf1=CCF_BINS; ccfmax(1)=X.SCNA_NA_CCF(i);
    ALT_DNA(1,:) = ccf1*purity*(X.SCNA_NA_CCF(i))*(X.SCNA_NA(i)>0);
    % 2) mutation on B allele NB, m=1
    ccf1=CCF_BINS; ccfmax(2)=X.SCNA_NB_CCF(i);
    ALT_DNA(2,:) = ccf1*purity*(X.SCNA_NB_CCF(i))*(X.SCNA_NB(i)>0);
    % 3) mutation on A allele NCNA, m=1
    ccf1=CCF_BINS; ccfmax(3)=1-X.SCNA_NA_CCF(i);
    ALT_DNA(3,:) = ccf1*purity*(X.SCNA_NCNA(i)*(1-X.SCNA_NA_CCF(i)));
    % 4) mutation on B allele NCNB, m=1
    ccf1=CCF_BINS; ccfmax(4)=1-X.SCNA_NB_CCF(i);
    ALT_DNA(4,:) = ccf1*purity*(X.SCNA_NCNB(i)*(1-X.SCNA_NB_CCF(i)));
    % 5) mutation on A allele NCNA and NA m=NA
    ccf1=CCF_BINS; ccfmax(5)=1;
    if (X.SCNA_NA(i)>0)
        ALT_DNA(5,:) = ccf1*purity*( (X.SCNA_NA(i)*X.SCNA_NA_CCF(i)) + (X.SCNA_NCNA(i)*(1-X.SCNA_NA_CCF(i))) );
    end
    % 6) mutation on B allele NCNB and NB m=NB
    ccf1=CCF_BINS; ccfmax(6)=1;
    if (X.SCNA_NB(i)>0)
        ALT_DNA(6,:) = ccf1*purity*( (X.SCNA_NB(i)*X.SCNA_NB_CCF(i)) + (X.SCNA_NCNB(i)*(1-X.SCNA_NB_CCF(i))) );
    end
    like1=0;
    model1=0;
    af1=0;
    P1=0*ccf1;
    for m1=1:6
        AF_DNA(m1,:)=ALT_DNA(m1,:)./totalDNA;
        AF_DNA(m1,AF_DNA(m1,:)>1)=0;
        % likelihood for ccfs + extended ccfs to allow for binominal AF fluctuation above max CCF
        ALT_LIKE_FULL(m1,:)=binopdf(X.t_alt_count(i),X.t_alt_count(i)+X.t_ref_count(i),AF_DNA(m1,:));
        % normalize
        %ALT_LIKE_FULL(m1,:)=ALT_LIKE_FULL(m1,:)/sum(ALT_LIKE_FULL(m1,:));
        % stdev
        [mu1,sigma1]=stats_mu_sigma(CCF_BINS,ALT_LIKE_FULL(m1,:) );
        % weight
        W=exp(-0.5*(CCF_BINS-ccfmax(m1,:)).^2./sigma1^2);
        W=W/max(W);
        W(CCF_BINS<=ccfmax(m1))=1;
        % restrict to possible ccf range
        %ALT_LIKE(m1,:)=ALT_LIKE_FULL(m1,:).*(CCF_BINS<=ccfmax(m1));
        % prob CCF conditional on CCF < ccfmax - attenutate ALT_LIKE
        P=ALT_LIKE_FULL(m1,:).*W;
        PX=P(ccf1>ccfmax(m1));
        kend=find(ccf1<=ccfmax(m1)); kend=kend(end);
        if (sum(PX)>0)
            P(kend)=P(kend)+sum(PX);
            P((kend+1):end)=0;
        end
        P=P.*CCF_PRIOR;
        %P=P/sum(P);
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
    X.ccf_ALT_model(i,1)=model1;
    X.ccf_ALT_model_AF(i,1)=af1;
    X.ccf_ALT_model_OVER(i)=pccfx/sum(pccf);
    
    if isnan(sum(pccf))
        strcat(X.Chromosome(i),':',num2str(X.Start_position(i)),'@',X.Tumor_Sample_Barcode(i),'@',X.Hugo_Symbol(i),'@',X.Variant_Classification(i))
        %pause()
        continue;
    end

    X.pccf(i,:)=pccf(CCF_BINS<=1)/sum(pccf(CCF_BINS<=1));
    
    %X.ccf_hat(i,1)=sum(ccf.*X.pccf(i,:))./sum(X.pccf(i,:));

    nb=length(AF_BINS);
    cccf=cumsum(X.pccf(i,:));
    ilow=sum(cccf<0.16);
    ihigh=sum(cccf<0.84)+1;
    if (ilow<1), ilow=1;end
    if (ihigh<=ilow), ihigh=ilow+1;end
    if (ihigh>nb), ihigh=nb;end
    X.ccf_CI_low(i,1)=CCF_BINS(ilow);
    X.ccf_CI_high(i,1)=CCF_BINS(ihigh);

    ilow=sum(cccf<0.025);
    ihigh=sum(cccf<0.975)+1;
    imedian=sum(cccf<0.5);
    if (imedian<1), imedian=1; end
    [~,imax]=max(pccf);
    if (ilow<1), ilow=1;end
    if (ihigh<=ilow), ihigh=ilow+1;end
    if (ihigh>nb), ihigh=nb;end
    X.ccf_CI95_low(i,1)=CCF_BINS(ilow);
    X.ccf_CI95_high(i,1)=CCF_BINS(ihigh);
    X.ccf_mode(i,1)=CCF_BINS(imax(1));
    X.ccf_mean(i,1)=sum(AF_BINS.*X.pccf(i,:))./sum(X.pccf(i,:));
    X.ccf_median(i,1)=CCF_BINS(imedian(1));
    X.multiplicity(i,1)=1;
    if model1==5
        X.multiplicity(i,1)=X.SCNA_NA(i);
    end
    if model1==6
        X.multiplicity(i,1)=X.SCNA_NB(i);
    end
    
    X.ccf_hat(i,:) = X.ccf_median(i,:);
    % set clonal flags
    X.clonal(i,1) = X.ccf_hat(i)>=0.9;
    
    %fprintf(1,' %d:%d\n',i,X.clonal(i,1))
end

%keyboard

if PAR.CONVERT_PCCF_TO_FIELDS
    for j=1:length(X.pccf(1,:))
        %vlab=sprintf('pccf_%d',round(100*ccf(j)));
        %X.(vlab)=X.pccf(:,j);
        vlab=sprintf('pccf_0z%d',round(100*ccf(j)));
        if (ccf(j)==0)
            vlab='pccf_0';
        end
        if (ccf(j)==1)
            vlab='pccf_1';
        end
        X.(vlab)=X.pccf(:,j);    
    end
end

if ~PAR.KEEP_PCCF_VECTOR
    ff=fieldnames(X);
    fx=ff(ismember(ff,{'paf','pccf','paft'}));
    if ~isempty(fx) 
        X=rmfield(X,fx);
    end
end

X1=X;
if ~isempty(PAR.OUTPUT_MAF)
    printStruct(X1,-1,PAR.OUTPUT_MAF);
end

if PAR.ADD_ABSOLUTE_FIELDS
    XA=v3toABSOLUTE_maf(X);
    if ~isempty(PAR.OUTPUT_MAF)
        ABS_MAF=regexprep(PAR.OUTPUT_MAF,'maf$','ABS.maf');
        printStruct(XA,-1,ABS_MAF);
    end
end

end

%%
function [mu,sig]=stats_mu_sigma(x,p)
 p=p/sum(p);
 mu=sum(x.*p);
 sig=sqrt( sum(x.^2.*p) - mu.^2);
end

%%
function [XA]=v3toABSOLUTE_maf(X3)
% add ABSOLUTE fields to V3 maf ... for phylogicNDT
XA=X3;
XA.observed_alt=X3.t_alt_count;
if isnumeric(XA.Chromosome)
    XA.Chromosome=num2str(XA.Chromosome);
    XA.Chromosome(find(ismember(XA.Chromosome,'23')))='X';
    XA.Chromosome(find(ismember(XA.Chromosome,'24')))='Y';
    XA.Chromosome(find(ismember(XA.Chromosome,'25')))='MT';
end
XA.male_X=ismember(upper(XA.Chromosome),'X')&(XA.SCNA_NCNB==0);
XA.A1_Z_ix=NaN*XA.SCNA_NCNA;  % probably SCNA cluster on A
XA.A2_Z_ix=NaN*XA.SCNA_NCNA;  % probably SCNA cluster on B
XA.NumberOfTimesCodonChangeIsInCOSMIC=XA.COSMIC_n_overlapping_mutations.*(cellfun(@length,XA.Protein_Change)>0);
XA.q_hat=XA.SCNA_NA+XA.SCNA_NB;
XA.HS_q_hat_1=XA.SCNA_NA;
XA.HS_q_hat_2=XA.SCNA_NB;
XA.SC_Z_ix=NaN*XA.SCNA_NCNA; % subclonal SCNA?
XA.C_Z_ix=NaN*XA.SCNA_NCNA; % clonal SCNA?
XA.total_qc=NaN*XA.SCNA_NCNA; % ?
XA.total_qs=NaN*XA.SCNA_NCNA; % ?
XA.SC_Aq_d=NaN*XA.SCNA_NCNA; % ?
XA.SC_Aq_a=NaN*XA.SCNA_NCNA; % ?
XA.C_Aq=NaN*XA.SCNA_NCNA; % ?
XA.clonal_scna_mut_ix=(XA.ccf_ALT_model==1).*(XA.SCNA_NA_CCF>0.9);
XA.Pr_somatic_clonal=sum(X3.pccf(:,90:end),2)
XA.Pr_germline=0*XA.Start_position;
XA.Pr_subclonal=1-sum(X3.pccf(:,90:end),2)
XA.Pr_subclonal_wt0=NaN*XA.SCNA_NCNA; % ?
XA.Pr_wt0=NaN*XA.SCNA_NCNA; % ?
XA.Pr_ge2=NaN*XA.SCNA_NCNA; % ?
XA.Pr_GL_som_HZ_alt=NaN*XA.SCNA_NCNA; % ?
XA.Pr_GL_som_HZ_ref=NaN*XA.SCNA_NCNA; % ?
XA.Pr_cryptic_SCNA=NaN*XA.SCNA_NCNA; % ?
XA.modal_q_s=XA.multiplicity; % ?
XA.LL=NaN*XA.SCNA_NCNA; % ?
XA.subclonal_Z_ix=XA.ccf_hat<0.9;
XA.Pr_subclonal_wt0=NaN*XA.SCNA_NCNA; % ?
XA.clonal_Z_ix=XA.ccf_hat>=0.9;
XA.wt0_Z_ix=NaN*XA.SCNA_NCNA; % ?
XA.clonal_het_Z_ix=1+0*XA.Start_position;
XA.ge2_Z_ix=NaN*XA.SCNA_NCNA; % ?
XA.homozygous_Z_ix=0*XA.Start_position;
XA.H1=NaN*XA.SCNA_NCNA; % ?
XA.H2=NaN*XA.SCNA_NCNA; % ?
XA.H3=NaN*XA.SCNA_NCNA; % ?
XA.H4=NaN*XA.SCNA_NCNA; % ?
XA.detection_power=NaN*XA.SCNA_NCNA; % ?
XA.detection_power_for_single_read=NaN*XA.SCNA_NCNA; % ?
XA.SSNV_skew=1+0*XA.SCNA_NCNA; % ?
ccf=0:0.01:1;
for j=1:length(X3.pccf(1,:))
    vlab=sprintf('pccf_0Z%d',round(100*ccf(j)));
    if (ccf(j)==0)
        vlab='pccf_0';
    end
    if (ccf(j)==1)
        vlab='pccf_1';
    end
    XA.(vlab)=X3.pccf(:,j);
end
end

%%
function test0
clear
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/Mutations/test')
X=load('../REBC_primary_pair_393.consensus_fix_NALT01.14Feb2020.maf.mat')
A=load_tsv('../../CN/CCF/REBC_primary_pair_393.alleliccapseg_pp_ccf_fit_v3_collapsed.tsv')
P=load_tsv('../../FC/REBC.pairs.15Mar2020.tsv')
PT=load_tsv('../../FC/REBC.participants.23Feb2020.tsv')

pid=unique(sort(A.sample))
P=trimStruct(P,ismember(P.entity_pair_id,pid))
PT=trimStruct(PT,ismember(PT.entity_participant_id,P.participant))

XV3=[]
for p=385:P.N
    pid1=pid{p};
    [num2str(p) '. ' pid1]
    A1=trimStruct(A,ismember(A.sample,pid1));
    P1=trimStruct(P,ismember(P.entity_pair_id,pid1));
    PT1=trimStruct(PT,ismember(PT.entity_participant_id,P1.participant));
    purity=A1.purity(1);
    ploidy=A1.ploidy(1);
    gender=PT1.Gender;
    X1=trimStruct(X,ismember(X.pID,P1.participant{1}));
    [ X3 ] = MAF_AC_PP_CCF_fit_v3(X1,A1,purity,ploidy,gender);
    if isempty(XV3)
        XV3=X3;
    else
        XV3=mergeStruct(XV3,X3);
    end
end    

save(['REBC_primary_pair_393.consensus_fix_NALT01.14Feb2020.CCF_fit_v3.' TODAY '.maf.mat'],'-struct','XV3')

end

%%
function test2
clear
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/Mutations/test')
X=load('../REBC_primary_pair_393.consensus_fix_NALT01.14Feb2020.maf.mat')
X3=load('REBC_primary_pair_393.consensus_fix_NALT01.14Feb2020.CCF_fit_v3.29Mar2020.maf.mat')
XA=load('../REBC_primary_pair_393.absolute_segforcecall_fix.16Feb2020.maf.mat')

X3.key=regexprep(strcat(X3.Tumor_Sample_Barcode,':',X3.Hugo_Symbol,':',X3.Chromosome,':',cellstr(num2str(X3.Start_position))), ' ','');
XA.key=regexprep(strcat(XA.Tumor_Sample_Barcode,':',XA.Hugo_Symbol,':',XA.Chromosome,':',cellstr(num2str(XA.Start_position))), ' ','');
[i3 mA]=ismember(X3.key,XA.key);
tab(i3)
X3.ABS_ccf_hat=NaN*X3.ccf_hat;
X3.ABS_ccf_hat(i3)=XA.ccf_hat(mA(i3));
X3.ABS_modal_q_s=NaN*X3.ccf_hat;
X3.ABS_modal_q_s(i3)=XA.modal_q_s(mA(i3));

hist2(X3.ABS_ccf_hat,X3.ccf_hat,100)
myhist(X3.ABS_ccf_hat-X3.ccf_hat,100,'log')
klow=find((X3.ABS_ccf_hat-X3.ccf_hat)<-0.15)
X3.Tumor_Sample_Barcode(klow)
q=trimStruct(X3,klow); q=rmfield(q,{'N','pccf','CCF_low','CCF_high'}); q=struct2table(q)
% AF8W chr 13 V3 ok - subclonal normal A allele 
khi=find((X3.ABS_ccf_hat-X3.ccf_hat)>0.2)
X3.Tumor_Sample_Barcode(khi)
q=trimStruct(X3,khi); q=rmfield(q,{'N','pccf','CCF_low','CCF_high'}); q=struct2table(q)
% AF8W chr 13 V3 ok - subclonal normal A allele 


end
function test
%% REBC 

clear
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/Mutations/test')
pid1='REBC-ACAB-TP-NB'
%pid1='REBC-ACAE-TP-NT'
X1=load_tsv([pid1 '_consensus_fix_NALT01.maf'])
A=load_tsv('../../CN/CCF/REBC_primary_pair_393.alleliccapseg_pp_ccf_fit_v3_collapsed.tsv')
P=load_tsv('../../FC/REBC.pairs.15Mar2020.tsv')
PT=load_tsv('../../FC/REBC.participants.23Feb2020.tsv')
A1=trimStruct(A,ismember(A.sample,pid1))
P1=trimStruct(P,ismember(P.entity_pair_id,pid1))
PT1=trimStruct(PT,ismember(PT.entity_participant_id,P1.participant))
purity=A1.purity(1)
ploidy=A1.ploidy(1)
gender=PT1.Gender;
[ X3 ] = MAF_AC_PP_CCF_fit_v3(X1,A1,purity,ploidy,gender)
[A1.CCF_hat,k]=max([A1.NA_CCF A1.NB_CCF],[],2)
A1.CCF_Low=A1.CCF_hat-0.1;
A1.CCF_High=1+0*A1.CCF_hat;
A1.NCN=2+0*A1.CCF_hat;
[ X2A ] = MAF_AC_PP_CCF_fit(X1,A1,purity,ploidy,gender)

clf
plot(X2A.CCF_hat,X3.ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')

[q,k]=min(X3.ccf_hat)
% 395


%NA, NB, CCF_hat

myhist(X3.multiplicity,1:10)
myhist(X3.ccf_hat,100)
myhist(X3.ccf_median,100)
myhist(X3.ccf_mode,100)

X2=load_tsv([pid1 '.CCF.maf'])
XA=load_tsv([pid1 '_ABS_MAF.txt'])
X3.key=regexprep(strcat(X3.Hugo_Symbol,':',X3.Chromosome,':',cellstr(num2str(X3.Start_position))), ' ','');
X2.key=regexprep(strcat(X2.Hugo_Symbol,':',X2.Chromosome,':',cellstr(num2str(X2.Start_position))), ' ','');
XA.key=regexprep(strcat(XA.Hugo_Symbol,':',cellstr(num2str(XA.Chromosome)),':',cellstr(num2str(XA.Start_position))), ' ','');
[i3 mA]=ismember(X3.key,XA.key);
tab(i3)
X3.ABS_ccf_hat=NaN*X2.CCF_hat;
X3.ABS_ccf_hat(i3)=XA.ccf_hat(mA(i3));
X3.ABS_modal_q_s=NaN*X2.CCF_hat;
X3.ABS_modal_q_s(i3)=XA.modal_q_s(mA(i3));
[i3 mA]=ismember(X3.key,X2.key);
tab(i3)
X3.V2_ccf_hat=NaN*X2.CCF_hat;
X3.V2_ccf_hat(i3)=X2.CCF_hat(mA(i3));
clf
plot(X3.ccf_hat,X3.ABS_ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--'); xlabel('V3 CCF'); ylabel('ABS CCF')
plot(X3.ccf_hat,X3.V2_ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--'); xlabel('V3 CCF'); ylabel('ABS CCF')
clf
hist2(X3.ccf_hat,X3.ABS_ccf_hat,100); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--'); xlabel('V3 CCF'); ylabel('ABS CCF')

clf
myhist(X3.ABS_ccf_hat,100)
myhist(X3.ccf_hat,100)

imagesc(X3.pccf)

X3.ccf_delta=X3.ccf_hat-X3.ABS_ccf_hat;
plot(X3.ccf_delta)
ff=fieldnames(XA)
ff=ff(strfindk(ff,'pccf_'))
z=[]
for j=1:length(ff)
    q=X2.(ff{j});
    if ~isnumeric(q)
        q=str2double(q);
    end
    z(:,j)=q';
end
clf
k=find(X3.ccf_delta<-0.162)
plot(X3.pccf(k,:))
k2=find(ismember(X2.key,X3.key(k)))
hold on
plot(z(k2,:))
hold off



plot(X3.ABS_ccf_hat,X3.ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')
plot(X3.ABS_ccf_hat,X3.ccf_median,'o')
plot(X3.ABS_ccf_hat,X3.ccf_mode,'o')
plot(X2A.CCF_hat,X3.ABS_ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')
plot(X2A.CCF_hat,X3.V2_ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')


%% check outliers in ACAB
plot(X3.ABS_ccf_hat,X3.ccf_hat,'o');line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--'); ylabel('V3 CCF'); xlabel('ABS CCF')
X3.ccf_delta=X3.ccf_hat-X3.ABS_ccf_hat;
myhist(X3.ccf_delta,100,'log')
tabulate(X3.Chromosome)
i=X3.ccf_delta<-0.6;
tabulate(X3.Chromosome(i))
clf
plot(X3.ccf_delta)
hold on;
plot(X3.SCNA_NA)
plot(X3.SCNA_NB)
plot(X3.SCNA_NA_CCF)
plot(X3.SCNA_NB_CCF)
plot(X3.i_tumor_f*10)
plot(X3.multiplicity+0.1,':')
plot(X3.ccf_ALT_model,':')

hold off


(X3.NA(i))



%%
hist2(X3.ABS_modal_q_s,X3.multiplicity,1:7)

[i3 m2]=ismember(X3.key,X2.key);
tab(i3)
X3.V2_ccf_hat=NaN*X2.CCF_hat;
X3.V2_ccf_hat(i3)=X2.CCF_hat(m2(i3));
X3.V2_ccf_median=NaN*X2.CCF_hat;
X3.V2_ccf_median(i3)=X2.CCF_median(m2(i3));
clf
plot(X3.V2_ccf_hat,X3.ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')
plot(X3.V2_ccf_median,X3.ccf_median,'o')

clf
plot(X3.ABS_ccf_hat,X3.V2_ccf_hat,'o'); line([0 1], [0 1],'color',0.5*[1 1 1],'linestyle','--')

hist2(X3.ABS_modal_q_s,X3.multiplicity,1:7)

hist2(X3.ABS_modal_q_s,X3.multiplicity,1:7)

myhist(XA.HS_q_hat_1,20)

% 
% 
% 
% P.sample=cellfun(@(x) x(1),regexp(P.pair_id,'-TP-','split'))
% P.sample=regexprep(P.sample,'-nullpair$','')
% [i m]=ismember(AC.pair_id,P.pair_id)
% AC.sample(i,1)=P.sample(m(i))
% a=tab(AC.sample)
% s=a.x{a.n==max(a.n)};
% s=a.x{a.n==min(a.n)};
% 
% XA=load_tsv('/xchip/cga_home/stewart/Cancer/DLBCL/ABSOLUTE/DLBCL_Paired_No_Artifact/ABSOLUTE_results/DLBCL/reviewed/SEG_MAF/DLBCL-LS1899-TP-NB-SM-1Q67X-SM-1Q65N_ABS_MAF.txt','trim')
% TA=load_tsv('/xchip/cga_home/stewart/Cancer/DLBCL/ABSOLUTE/DLBCL_Paired_No_Artifact/ABSOLUTE_results/DLBCL/reviewed/SEG_MAF/DLBCL-LS1899-TP-NB-SM-1Q67X-SM-1Q65N.segtab.txt')
% % pccf struct
% ff=fieldnames(XA)
% k=find(cellfun(@length,regexp(ff,'^x\d+','match'))>0)
% ff=ff(k);
% n0=length(ff)
% XA.ccf_ABS=0:(1/(n0-1)):1;
% for j=1:n0
%     XA.pccf_ABS(:,j)=XA.(ff{j});
% end
% s=upper(regexprep(XA.sample{1},'-Tumor',''))
% AC1=trimStruct(AC,ismember(AC.sample,s))
% purity1=P.purity(ismember(P.sample,s))
% ploidy1=P.ploidy(ismember(P.sample,s))
% gender1=G.sex{ismember(G.sample,s)}
% [i m]=ismember(P.sample,s)
% X1=trimStruct(X,ismember(X.Tumor_Sample_Barcode,P.case_sample(i)))
% 
% T1=AllelicCapSeg_PP_CCF_fit(AC1,purity1,ploidy1,gender1)
% X1= MAF_AC_PP_CCF_fit(X1, T1 ,purity1, ploidy1, gender1)
% 
% X1.id=regexprep(strcat(X1.Chromosome,':',cellstr(num2str(X1.Start_position))), ' ','');
% XA.id=regexprep(strcat(num2chrom(XA.Chromosome),':',cellstr(num2str(XA.Start_position))), ' ','');
% 
% [i m]=ismember(X1.id,XA.id);
% X1m=trimStruct(X1,find(i))
% XAm=trimStruct(XA,m(i))
% plot(XAm.ccf_hat,X1m.CCF_hat,'+')
% hist2(XAm.q_hat,X1m.SCNA_q_hat,100)
% 
% X1M= MAF_AC_PP_CCF_fit(X1m, T1 ,purity1, ploidy1, gender1)
% 
% T1.id=regexprep(strcat(num2chrom(T1.Chromosome),':',cellstr(num2str(T1.Start_bp))), ' ','');
% TA.id=regexprep(strcat(num2chrom(TA.Chromosome),':',cellstr(num2str(TA.Start_bp))), ' ','');
% 
% TA.cancer_cell_frac=max([str2double(TA.cancer_cell_frac_a1) str2double(TA.cancer_cell_frac_a2)],[],2)
% [i m]=ismember(T1.id,TA.id);
% T1m=trimStruct(T1,find(i))
% TAm=trimStruct(TA,m(i))
% plot(TAm.cancer_cell_frac,T1m.CCF_hat,'+')
% 
% CCF_hat
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

