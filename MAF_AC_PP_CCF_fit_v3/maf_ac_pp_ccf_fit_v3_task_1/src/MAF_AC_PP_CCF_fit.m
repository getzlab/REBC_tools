function [ X1 ] = MAF_AC_PP_CCF_fit(X,T,purity,ploidy,gender,PAR)
%  [ X1 ] = MAF_AC_PP_CCF_fit(X, T ,purity,ploidy,gender,PAR)
%
%  MAF_AC_PP_CCF_fit adds estimated mutation CCF and local tumor SCNA allele copy number 
%  [NA,NB] and  CCF's to AC based on puity and ploidy (from ABSOLUTE or other method).
%  This function makes estimates for all segments, even X and Y when gender
%  info is available from PAR. 
%  MAF_AC_PP_CCF_fit operates on one sample at a time (one purity, ploidy, gender) 
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
%     T: ALLELIC_CAPSEG modified struct w/ fitted (NA,NB,CCF_hat) fields: 
%       tau, sigma_tau, mu_major, sigma_major, mu_minor,sigma_minor, NA, NB, CCF_hat
%     purity: tumor sample fraction of cancer cells (0 to 1)
%     ploidy: cancer cell average CN over genome (0 to  inf)
%     gender: 'MALE' or 'FEMALE' ('M' or 'F' in any case combination) 
%     PAR fields:
%       MAXCN: maximum tumor SCNA copy number (10)
%       SMOOTH_SPAN: taps in chisquare smoothing function (7)
%       CCF_BINS: bins over CCF range (0-1)
%       SIGMA_THRESHOLD: signal:noise threshold for SCNA (2)
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
    PAR=[]
    PAR.DEFAULT=1;
end
if ~isfield(PAR,'CCF_BINS')
    PAR.CCF_BINS=0:0.01:1; 
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
    PAR.ASSUME_NORMAL_X_WHEN_MISSING=0; 
end
XCN=2; % X normal CN 
if ismember(upper(gender),{'M','MALE'})
    XCN=1;
end


X.SCNA_NA=NaN*X.Start_position;
X.SCNA_NB=NaN*X.Start_position;
X.SCNA_q_hat=NaN*X.Start_position;
X.SCNA_CCF_hat=NaN*X.Start_position;
X.SCNA_CCF_Low=NaN*X.Start_position;
X.SCNA_CCF_High=NaN*X.Start_position;
X.SCNA_tau=NaN*X.Start_position;
X.SCNA_NCN=NaN*X.Start_position;
X.IS_SCNA=NaN*X.Start_position;

X.purity=purity*ones(size(X.Start_position));
X.ploidy=ploidy*ones(size(X.Start_position));

X.dna_fraction_in_tumor=NaN*X.Start_position;
X.CCF_hat=NaN*X.Start_position;
X.CCF_Low=NaN*X.Start_position;
X.CCF_High=NaN*X.Start_position;

X.CCF_CI95_low=NaN*X.Start_position;
X.CCF_CI95_high=NaN*X.Start_position;
X.CCF_mode=NaN*X.Start_position;
X.CCF_mean=NaN*X.Start_position;
X.CCF_median=NaN*X.Start_position;
X.i_tumor_f_corrected=NaN*X.Start_position;
X.pccf=NaN(length(X.Start_position),length(PAR.CCF_BINS));
X.clonal = NaN*X.Start_position;

% af binning
ccf=PAR.CCF_BINS;
daf=median(diff(PAR.CCF_BINS));
% alt allele fraction
if ~isfield(X,'i_tumor_f')
    X.i_tumor_f=X.t_alt_count./(X.t_alt_count+X.t_ref_count);
end

if isnumeric(X.Chromosome)
    X.Chromosome=num2chrom(X.Chromosome);
end
    
N=length(X.Start_position)
for i=1:N
    C=X.Chromosome(i);
    if ~isnumeric(X.Chromosome)
        C=chrom2num(C);
    end
    
    k=find(ismember(T.Chromosome,C)&(T.Start_bp<=X.Start_position(i))&(T.End_bp>=X.Start_position(i)));
    if numel(k)==0 
        % no SCNA segment
        if (C>22)&PAR.ASSUME_NORMAL_X_WHEN_MISSING
            % X or Y
            NCN=0;
            if C==23
                NCN=XCN;
            elseif C==24
                NCN=1-XCN;
            end
            if (NCN<1), continue; end  % no Y
            k=length(T.NA)+1;
            T.NCN(k)=NCN;
            T.NA(k)=0;
            T.NB(k)=1;
            T.CCF_hat(k)=0;
            T.CCF_Low(k)=0;
            T.CCF_High(k)=1;
            T.tau(k)=NaN;
            T.IS_SCNA(k)=NaN;
        else
            continue;
        end
    end
    % if two segments overlap (problem) pick 1 
    if numel(k)>1, k=k(1); end
    % assign SCNA state at mutation site
    X.SCNA_NA(i)=T.NA(k);
    X.SCNA_NB(i)=T.NB(k);
    X.SCNA_CCF_hat(i)=T.CCF_hat(k);
    X.SCNA_CCF_Low(i)=T.CCF_Low(k);
    X.SCNA_CCF_High(i)=T.CCF_High(k);
    X.SCNA_tau(i)=T.tau(k);
    X.IS_SCNA(i)=T.IS_SCNA(k);
    X.SCNA_NCN(i)=T.NCN(k);    
    X.SCNA_q_hat(i)=X.SCNA_NA(i)+X.SCNA_NB(i);
    tumorDNA=X.SCNA_q_hat(i).*X.SCNA_CCF_hat(i)+X.SCNA_NCN(i).*(1-X.SCNA_CCF_hat(i));
    X.dna_fraction_in_tumor(i)=(X.purity(i).*tumorDNA)./(X.purity(i).*tumorDNA+X.SCNA_NCN(i).*(1-X.purity(i)));
    
    % allele fraction pdf
    X.paf(i,:)=daf*betapdf(PAR.CCF_BINS,X.t_alt_count(i)+1,X.t_ref_count(i)+1);
    % x1.paf(j,:) not always normalized (t_alt_count = 0 problem)
    if abs(sum(X.paf(i,:))-1)>0.1
        'paf',sum(X.paf(i,:)) % keyboard;
    end
    X.paf(i,:)=X.paf(i,:)/sum(X.paf(i,:));
    X.pccf(i,:)=NaN*X.paf(i,:);
    X.paft(i,:)=NaN*X.paf(i,:);
    % skip events that ABSOLUTE skipped
    if isnan(X.dna_fraction_in_tumor(i)), continue; end
    % all possible ref count in tumor from 0 to t_ref_count
    nreft=0:X.t_ref_count(i);
    % likelihood for each possible tumor ref count
    wreft=binopdf(nreft+X.t_alt_count(i),X.t_ref_count(i)+X.t_alt_count(i),X.dna_fraction_in_tumor(i));
    % normalize wreft within possible nreft counts
    wreft=wreft/sum(wreft);
    % init tumor specific allele fraction dist
    paft=0*PAR.CCF_BINS; paftv=[];
    % loop over each tumor ref count
    for n=1:length(nreft)
        paft1=daf*wreft(n)*betapdf(PAR.CCF_BINS,X.t_alt_count(i)+1,nreft(n)+1);
        paftv=[paftv; paft1];
        paft=paft+paft1;
    end
    % store tumor allele fracton pdf to maf
    if (abs(sum(paft)-1)>0.1)&(purity>0.8)
        'paft',sum(paft) % keyboard;
    end
    X.paft(i,:)=paft/sum(paft);
    % possible multiplicities depending on when the mutation happened relative to the SCNA
    m=unique(sort([1 X.SCNA_NA(i) X.SCNA_NB(i)]));
    m=m(m>0);
    % half weight for m=1 (after SNCA), rest for m>1 for each m (before SCNA)
    wm=[];
    % prior SCNA ccf to bias mutation to mult=1 for high ccf SCNAs (SCNA before mutation)
    wm(m==1,1)=1e-10+X.SCNA_CCF_hat(i);
    if any(m>1)
        wm(m>1,1)=(1-X.SCNA_CCF_hat(i))/sum(m>1);
    end
    wm=wm/sum(wm);
    % correct for SCNA ccf
    % mutation can appear on SCNA_ccf_hat or copy 1 with 1-SCNA_ccf_hat
    wm(1)=1-X.SCNA_CCF_hat(i)*(1-wm(1));
    wm(2:end)=X.SCNA_CCF_hat(i)*wm(2:end);
    wm/sum(wm);
    %keyboard
    %ccf=0:daf:1;
    ccf=PAR.CCF_BINS; %(daf/2):daf:1;
    pccf=0*ccf;
    pccfv=[];
    m=1;
    for m=1:length(m)
        % scale ccf down by n
        ccfb=(X.SCNA_NCN(i)*(1-X.SCNA_CCF_hat(i))+X.SCNA_CCF_hat(i)*X.SCNA_q_hat(i))*(ccf/m(m));
        % project aft into ccfb regions with ccf bins
        if (PAR.SMOOTH>1)
            ccfx=0:(1/PAR.SMOOTH):1;
            px=interp1(ccf,X.paft(i,:),ccfx);
            ccfb=(X.SCNA_NCN(i)*(1-X.SCNA_CCF_hat(i))+X.SCNA_CCF_hat(i)*X.SCNA_q_hat(i))*(ccfx/m(m));
            p1=histw(px,ccfb,ccf);
            p1=p1*length(PAR.CCF_BINS)/PAR.SMOOTH;
        else    
            p1=histw(X.paft(i,:),ccfb,ccf);
        end
        %keyboard;
        % weight by wm
        pccf1=wm(m)*p1;
        pccf=pccf+pccf1;
        pccfv=[pccfv; pccf1];
    end
    pccfx=1-sum(pccf);
    %pccf=pccf/sum(pccf);
    X.pccf(i,:)=pccf;
    if (pccfx>0)
        X.pccf(i,end)=X.pccf(i,end)+1-sum(pccf);
    end
    if (X.pccf(i,end)<0), X.pccf(i,end)=0;end
    if isnan(sum(pccf))
        strcat(X.Chromosome(i),':',num2str(X.Start_position(i)),'@',X.Tumor_Sample_Barcode(i),'@',X.Hugo_Symbol(i),'@',X.Variant_Classification(i))
        %pause()
        continue;
    end
    nb=length(PAR.CCF_BINS);
    cccf=cumsum(X.pccf(i,:));
    ilow=sum(cccf<0.16);
    ihigh=sum(cccf<0.84)+1;
    if (ilow<1), ilow=1;end
    if (ihigh<=ilow), ihigh=ilow+1;end
    if (ihigh>nb), ihigh=nb;end
    X.CCF_CI_low(i,1)=ccf(ilow);
    X.CCF_CI_high(i,1)=ccf(ihigh);

    ilow=sum(cccf<0.025);
    ihigh=sum(cccf<0.975)+1;
    imedian=sum(cccf<0.5);
    if (imedian<1), imedian=1; end
    [pmax,imax]=max(pccf);
    if (ilow<1), ilow=1;end
    if (ihigh<=ilow), ihigh=ilow+1;end
    if (ihigh>nb), ihigh=nb;end
    X.CCF_CI95_low(i,1)=ccf(ilow);
    X.CCF_CI95_high(i,1)=ccf(ihigh);
    %x1.cancer_cell_frac(i,1)=ccf(imax);
    X.CCF_mode(i,1)=ccf(imax(1));
    X.CCF_mean(i,1)=sum(ccf.*X.pccf(i,:))./sum(X.pccf(i,:));
    X.CCF_median(i,1)=ccf(imedian(1));
    [q k]=max(X.paft(i,:));
    X.i_tumor_f_corrected(i,1)=PAR.CCF_BINS(k);
    
    %y=sum(reshape(X.pccf(i,:),10,100),1);
    %xb=mean(reshape(ccf,10,100),1);
    %y1 = interp1(xb,y,ccfA,'linear','extrap'); y1(y1<0)=0; y1=y1/sum(y1);
    %X.pCCF(i,:)=y1;
    %X.pCCF(i,:)=pccf;
    [ymax,imode]=max(pccf);
    X.CCF_mode(i,1)=ccf(imode);
    
    % set clonal flags
    X.clonal(i,1) = X.CCF_median(i)>=0.85;
    
    fprintf(1,' %d:%d',i,X.clonal(i,1))
end

X.CCF_hat = X.CCF_median;
%keyboard

if PAR.CONVERT_PCCF_TO_FIELDS
    for j=1:length(X.pccf(1,:))
        vlab=sprintf('pccf_%d',round(100*ccf(j)))
        X.(vlab)=X.pccf(:,j);
    end
end

if ~PAR.KEEP_PCCF_VECTOR
    ff=fieldnames(X)
    fx=ff(ismember(ff,{'paf','pccf','paft'}));
    if ~isempty(fx) 
        X=rmfield(X,fx);
    end
end

X1=X;

if length(PAR.OUTPUT_MAF)>0
    printStruct(X1,-1,PAR.OUTPUT_MAF);
end

end
function test
%% DLBCL 

clear
cd /Users/stewart/Projects/Cancer/DLBCL/maf/
X=load_maf('/Users/stewart/Projects/Cancer/DLBCL/maf/DLBCL_All.22-Mar-2015.maf','trim')
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

XA=load_tsv('/xchip/cga_home/stewart/Cancer/DLBCL/ABSOLUTE/DLBCL_Paired_No_Artifact/ABSOLUTE_results/DLBCL/reviewed/SEG_MAF/DLBCL-LS1899-TP-NB-SM-1Q67X-SM-1Q65N_ABS_MAF.txt','trim')
TA=load_tsv('/xchip/cga_home/stewart/Cancer/DLBCL/ABSOLUTE/DLBCL_Paired_No_Artifact/ABSOLUTE_results/DLBCL/reviewed/SEG_MAF/DLBCL-LS1899-TP-NB-SM-1Q67X-SM-1Q65N.segtab.txt')
% pccf struct
ff=fieldnames(XA)
k=find(cellfun(@length,regexp(ff,'^x\d+','match'))>0)
ff=ff(k);
n0=length(ff)
XA.ccf_ABS=0:(1/(n0-1)):1;
for j=1:n0
    XA.pccf_ABS(:,j)=XA.(ff{j});
end
s=upper(regexprep(XA.sample{1},'-Tumor',''))
AC1=trimStruct(AC,ismember(AC.sample,s))
purity1=P.purity(ismember(P.sample,s))
ploidy1=P.ploidy(ismember(P.sample,s))
gender1=G.sex{ismember(G.sample,s)}
[i m]=ismember(P.sample,s)
X1=trimStruct(X,ismember(X.Tumor_Sample_Barcode,P.case_sample(i)))

T1=AllelicCapSeg_PP_CCF_fit(AC1,purity1,ploidy1,gender1)
X1= MAF_AC_PP_CCF_fit(X1, T1 ,purity1, ploidy1, gender1)

X1.id=regexprep(strcat(X1.Chromosome,':',cellstr(num2str(X1.Start_position))), ' ','');
XA.id=regexprep(strcat(num2chrom(XA.Chromosome),':',cellstr(num2str(XA.Start_position))), ' ','');

[i m]=ismember(X1.id,XA.id);
X1m=trimStruct(X1,find(i))
XAm=trimStruct(XA,m(i))
plot(XAm.ccf_hat,X1m.CCF_hat,'+')
hist2(XAm.q_hat,X1m.SCNA_q_hat,100)

X1M= MAF_AC_PP_CCF_fit(X1m, T1 ,purity1, ploidy1, gender1)

T1.id=regexprep(strcat(num2chrom(T1.Chromosome),':',cellstr(num2str(T1.Start_bp))), ' ','');
TA.id=regexprep(strcat(num2chrom(TA.Chromosome),':',cellstr(num2str(TA.Start_bp))), ' ','');

TA.cancer_cell_frac=max([str2double(TA.cancer_cell_frac_a1) str2double(TA.cancer_cell_frac_a2)],[],2)
[i m]=ismember(T1.id,TA.id);
T1m=trimStruct(T1,find(i))
TAm=trimStruct(TA,m(i))
plot(TAm.cancer_cell_frac,T1m.CCF_hat,'+')

CCF_hat
end
