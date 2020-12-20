function [X]=alleleFraction_clonal_hets_Purity_estimate(f,DEPTHZ,plot_area,DAFB,DRIVERS)

if nargin<2
    DEPTHZ=2 % DEPTH outlier removal default number of sigmas way from normal fit peak to remove
end
if nargin<3
    plot_area='' % print plots to this area: file name = alleleFraction_clonal_hets.<Tumor_Sample_Barcode>.af.png
end
if nargin<4
    DAFB=1/200 % allele fraction bin size
end
if nargin<5
    DRIVERS={}
end

VC={'Missense_Mutation' 'Splice_Site' 'Nonsense_Mutation' 'Start_Codon_SNP' 'Nonstop_Mutation','De_novo_Start_InFrame','De_novo_Start_OutOfFrame'}

PLOT=length(plot_area)>0
if (PLOT)
    default_area=pwd();
    full_plot_area=[default_area '/' plot_area];
    mkdir(full_plot_area);
    
end

if ischar(f)
    M=load_table(f,char(9));
else
    M=f;
end
% chromosomes 1-22
if (iscell(M.Chromosome))
    a=chrom2num(M.Chromosome);
else
    a=M.Chromosome;
end
M=trimStruct(M,find(a<23));

M=trimStruct(M,strfindk(M.Variant_Type,'SNP'));
x=M.t_alt_count;
if (iscell(x))
    M.t_alt_count=str2num(char(M.t_alt_count));
    M.t_ref_count=str2num(char(M.t_ref_count));
end

% samples
sample=unique(M.Tumor_Sample_Barcode);


% allele fraction bins
afb=0:DAFB:1;

% output struct
X.sample=sample;
X.afb=repmat(afb,length(sample),1);
X.paf=0*X.afb;
X.afhist=0*X.afb;
X.afbeta=0*X.afb;
X.median_depth=NaN*zeros(length(sample),1);
X.Nmut0=NaN*zeros(length(sample),1);
X.afp=NaN*zeros(length(sample),1);
X.pur=NaN*zeros(length(sample),1);
X.Nmut=zeros(length(sample),1);
X.afCI=NaN*zeros(length(sample),2);
X.purCI=NaN*zeros(length(sample),2);
X.nClonalMut=NaN*zeros(length(sample),1);

for i=1:length(X.sample)
    i
    s1=X.sample{i}
    % mutations in this sample
    ks=strfindk(M.Tumor_Sample_Barcode,s1);
    if (isempty(ks)),continue; end

    X.Nmut0(i)=length(ks);
    
    X.median_depth(i)=nanmedian(M.t_alt_count(ks)+M.t_ref_count(ks));

    depth=M.t_alt_count(ks)+M.t_ref_count(ks);
    [mu sig]=normfit(depth);
   
    kout=((depth<(mu-DEPTHZ*sig))+(depth>(mu+DEPTHZ*sig)));
    ks=ks(~kout);
    X.Nmut(i)=length(ks);
    
    % allele fraction histogram (no smoothing)
    X.afhist(i,:)=hist(M.t_alt_count(ks)./(M.t_alt_count(ks)+M.t_ref_count(ks)),afb);

    % binomial smoothing
    n=0;
    for k1=ks'    
        n=n+1;
        %beta pdf for each mutation
        x2=betapdf(afb,M.t_alt_count(k1)+1,M.t_ref_count(k1)+1);
        % add normalized beta pdf for each mutation
        X.beta_pdfs(i,1).Hugo_Symbol(n,1)=M.Hugo_Symbol(k1);
        X.beta_pdfs(i,1).Chromosome(n,1)=M.Chromosome(k1);
        X.beta_pdfs(i,1).Start_position(n,1)=M.Start_position(k1);
        X.beta_pdfs(i,1).Variant_Classification(n,1)=M.Variant_Classification(k1);
        X.beta_pdfs(i,1).t_alt_count(n,1)=M.t_alt_count(k1);
        X.beta_pdfs(i,1).t_ref_count(n,1)=M.t_ref_count(k1);
        X.beta_pdfs(i,1).pdf(n,:)=(x2/sum(x2)+1e-40);
        [~,kp1]=max(X.beta_pdfs(i,1).pdf(n,:));
        X.beta_pdfs(i,1).pdf_max(n,1)=afb(kp1);
        X.beta_pdfs(i,1).DRIVERGENE(n,1)={''};
        if length(DRIVERS)>0
            kd=find(ismember(M.Hugo_Symbol(k1),DRIVERS).*ismember(M.Variant_Classification(k1),VC));
            X.beta_pdfs(i,1).DRIVER(n,1)=length(kd);
            if length(kd)>0
                X.beta_pdfs(i,1).DRIVERGENE{n,1}=[M.Hugo_Symbol{k1} ':' M.Protein_Change{k1}];
            end
        end
        X.paf(i,:)=X.paf(i,:)+(x2/sum(x2));
        
    end
    % plot(X.beta_pdfs(1).pdf','b-')
    % this sample's af dist
    x1=X.paf(i,:);
    dx2=-diff([0 0 x1],2)/DAFB;
    %plot(afb,[x1 ;dx2]); grid on; 
    
    % find het peak below 0.5 
    j=2:(length(x1)-1);
    kp0=find((x1(j-1)<=x1(j))&(x1(j+1)<x1(j)));   
    % second derivative peak
    kp=find((dx2(j-1)<=dx2(j))&(dx2(j+1)<dx2(j)));    
    % second derivative valley
    kv=find((dx2(j-1)>dx2(j))&(dx2(j+1)>dx2(j)));    
    % peak below AF 0.5
    kp(afb(kp)>=0.5)=[];
    if isempty(kp)
        continue;
    end
    [kp0(end) kp(end)];
    xp=x1(kp+1);
    ap=afb(kp+1);
    ap=ap(ap<=0.5);
    xp=xp(ap<=0.5);
    ap=ap(xp>0.0025);
    xp=xp(xp>0.0025);
    [ap,k]=max(ap);
    xp=xp(k)
    %ap=min([0.5 ap]);
    if (~isempty(k))
        X.afp(i)=ap;
        X.pur(i)=2*X.afp(i);    
        % CI on purity (assuming it's the het peak)
        laf=log(X.paf(i,:)/sum(X.paf(i,:))+1e-5);
        b0=find(afb==ap);
        inpeak=(X.beta_pdfs(i).pdf(:,b0)>(DAFB)); kclonal=find(inpeak);ksubclonal=find(~inpeak);
   
        %bL=kv(find(kv<b0)); bL=bL(end);
        %bH=kv(find(kv>b0)); bH=bH(1);
        %q=X.beta_pdfs(i,1).pdf_max;
        %k1=find( (q>=afb(bL)& (q<=afb(bH)) ));
        X.nClonalMut(i)=length(kclonal);
        if length(kclonal)>0
            if length(kclonal)>1
                logLike=sum(log10(X.beta_pdfs(i,1).pdf(kclonal,:)));
            else
                logLike=(log10(X.beta_pdfs(i,1).pdf(kclonal,:)));
            end    
            logLike=logLike-max(logLike);
            % revise peak for clonal mutations only
            [x,b0]=max(logLike)
            
            kL=find(logLike(1:(b0-1))<-1); kL=kL(end)-1;
            kH=b0+find(logLike((b0+1):end)<-1); kH=kH(1)+1;
        else
            kL=1;
            kH=length(afb);
        end
        X.afp(i)=afb(b0);
        X.pur(i)=2*X.afp(i);    
  
        if kL>=b0, kL=b0-1; end
        if kH<=b0, kH=b0+1; end
        X.afCI(i,1)=max([0 (afb(kL))]);
        X.afCI(i,2)=min([1 (afb(kH))]);
%         lafT=laf(b0)-0.5;
%         for b=b0:length(afb)
%             if (laf(b)<=(lafT))
%                 X.afCI(i,2)=min([1 (afb(b))]);
%                 break;
%             end
%         end
%         for b=b0:-1:1
%             if (laf(b)<=(lafT))
%                 X.afCI(i,1)=(afb(b));
%                 break;
%             end
%         end
    end
    if (X.afCI(i,2)>1), X.afCI(i,2)=1;,end
    X.purCI(i,:)=2*X.afCI(i,:);
    if (X.purCI(i,2)>1), X.purCI(i,2)=1;,end
    
    %inpeak=(X.beta_pdfs(i).pdf(:,b0)>(DAFB)); kclonal=find(inpeak);
    %pafj=prod(X.beta_pdfs(i).pdf(kclonal,:))
        
    if (PLOT)
       clf
       pause(0.001)
       plot(afb,exp(laf),'k-','linewidth',2);
       set(gcf,'Position', [10 500 900 400]);
       hold on
       plot(1,1,'--','color',0.05*[1 1 1]);
       plot(1,1,'-','color',0.5*[1 1 1]);
       plot(1,1,'-','color',0.75*[1 1 1]);
       plot(1,1,'r--');
       plot(1,1,'b--');
       %inpeak=(X.beta_pdfs(i).pdf(:,b0)>(DAFB/10)); kclonal=find(inpeak); ksubclonal=find(~inpeak);
       skip=1
       if length(ks)>3000
           skip=round(length(ks)/100);
          % kclonal=kclonal(1:skip:end);
           ksubclonal=ksubclonal(1:skip:end);
       end
       plot(afb,X.beta_pdfs(i).pdf(kclonal,:)'/sqrt(length(ks))+1e-5,'-','color',0.5*[1 1 1]);
       %pafc=sum(X.beta_pdfs(i).pdf(kclonal,:));
       %pafc=(pafc/max(pafc))*exp(laf(b0));
       %plot(afb,pafc+1e-5,'--','color',0.05*[1 1 1]);
       if length(kclonal)>1
           pafc=sum(log10(X.beta_pdfs(i).pdf(kclonal,:)));
           pafc=pafc-max(pafc);
           pafc=10.^pafc;
       else
           pafc=(X.beta_pdfs(i).pdf(kclonal,:));
       end
       pafc=(pafc/max(pafc))*exp(laf(b0));
       plot(afb,pafc+1e-5,'-','color',0.9*[1 0 0]);
      
       set(gca,'yscale','log')
       if length(ksubclonal)>0
           plot(afb,(X.beta_pdfs(i).pdf(ksubclonal,:)'/sqrt(length(ks))+1e-5),'-','color',0.75*[1 1 1]);
       end
       kdriver=find(X.beta_pdfs(i,1).DRIVER);
       txtdriver='';
       if length(kdriver)>0
           plot(afb,(X.beta_pdfs(i).pdf(kdriver,:)'/sqrt(length(ks))+1e-5),'--','color',0.75*[0 0 1]);
           txtdriver= X.beta_pdfs(i,1).DRIVERGENE{kdriver};
       end
            
       YL=exp([-11;-2.5]); ylim(YL); xlim([0 1])
       hL=line([1;1]*X.afp(i)',YL,'color','r','linestyle','--','linewidth',2);
       afLO=X.afCI(i,1);
       afHI=X.afCI(i,2);
       hL1=line([afLO;afHI],YL(2)*0.98*[1 1],'color','r','linestyle','-');
       hL1L=line(afLO*[1;1],YL(2)*0.98*[0.75; 1],'color','r','linestyle','-');
       hL1H=line(afHI*[1;1],YL(2)*0.98*[.75; 1],'color','r','linestyle','-');
       plot(afb,exp(laf),'k-','linewidth',2);
       
       %plot(afb,10.^(logLike)/100,'r-','linewidth',2);
       
       hold off
       text(0.96,0.95,s1,'units','normalized','fontsize',14,'horizontalalignment','right');
       skp=''
       if skip>1
           skp=[', skip=' num2str(skip)]
       end
       text(0.96,0.9,sprintf('n=%d, nclonal=%d, pur=%.2f [%.2f - %.2f] %s',X.Nmut(i),X.nClonalMut(i),X.pur(i),X.purCI(i,:),skp),'units','normalized','fontsize',14,'horizontalalignment','right');
       text(0.96,0.85,txtdriver,'units','normalized','fontsize',14,'horizontalalignment','right');
       if length(kdriver)>0
           h=legend({'cumulative pdf','cumulative clonal pdf','clonal het mutations','non-clonal het mutations','clonal het AF','driver'},'location','E');
       else
           h=legend({'cumulative pdf','cumulative clonal pdf','clonal het mutations','non-clonal het mutations','clonal het AF'},'location','E');    
       end

       %h=legend({'cumulative beta pdf','clonal mutations','non-clonal mutations','clonal het AF','driver'},'location','E');
       h.Position(2)=0.52;
       xlabel('Alt allele fraction','fontsize',16);
       ylabel('beta alt AF (normalized)','fontsize',16);
       grid on
       f=sprintf('%s/AF_purity_%s.%s',full_plot_area,s1,TODAY);
       print(gcf,'-dpng','-r70',[f '.png']);
       print(gcf,'-depsc',[ f '.eps']);

    end
end

   
    
end

