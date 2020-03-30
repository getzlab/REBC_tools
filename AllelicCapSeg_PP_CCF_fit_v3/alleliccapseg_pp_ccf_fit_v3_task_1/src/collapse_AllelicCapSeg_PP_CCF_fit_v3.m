function [x2]=collapse_AllelicCapSeg_PP_CCF_fit_v3(x1)

x1.x1=xhg19(x1.Chromosome,x1.Start_bp);
x1.x2=xhg19(x1.Chromosome,x1.End_bp);
[~,k]=sort(x1.x1);
x1=trimStruct(x1,k);
x2=trimStruct(x1,1);
i2=1;
for i1=2:x1.N
    if (x1.Chromosome(i1-1)==x1.Chromosome(i1))&(x1.cn_clusterA(i1-1)==x1.cn_clusterA(i1))&(x1.cn_clusterB(i1-1)==x1.cn_clusterB(i1))
        x2.End_bp(i2)=x1.End_bp(i1);
        n_probes2=x2.n_probes(i2);
        x2.n_probes(i2)=x1.n_probes(i1)+x2.n_probes(i2);
        len2=x2.length(i2);
        x2.length(i2)=x2.End_bp(i2)-x2.Start_bp(i2);
        x2.f(i2)=(x1.n_hets(i1)*x1.f(i1)+x2.n_hets(i2)*x2.f(i2))/(x1.n_hets(i1)+x2.n_hets(i2));
        x2.n_hets(i2)=x1.n_hets(i1)+x2.n_hets(i2);
        x2.tau(i2)=(x1.tau(i1)*x1.n_probes(i1)+x2.tau(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        % new sigma_tau should scale as standard error
        x2.sigma_tau(i2)=(x1.sigma_tau(i1)*x1.n_probes(i1)+x2.sigma_tau(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.mu_minor(i2)=(x1.mu_minor(i1)*x1.n_probes(i1)+x2.mu_minor(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.sigma_minor(i2)=(x1.sigma_minor(i1)*x1.n_probes(i1)+x2.sigma_minor(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.mu_major(i2)=(x1.mu_major(i1)*x1.n_probes(i1)+x2.mu_major(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.sigma_major(i2)=(x1.sigma_major(i1)*x1.n_probes(i1)+x2.sigma_major(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        
        x2.cn_minor(i2)=(x1.cn_minor(i1)*x1.n_probes(i1)+x2.cn_minor(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.cn_sigma_minor(i2)=(x1.cn_sigma_minor(i1)*x1.n_probes(i1)+x2.cn_sigma_minor(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.cn_major(i2)=(x1.cn_major(i1)*x1.n_probes(i1)+x2.cn_major(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.cn_sigma_major(i2)=(x1.cn_sigma_major(i1)*x1.n_probes(i1)+x2.cn_sigma_major(i2)*n_probes2)/(x1.n_probes(i1)+n_probes2);
        x2.x2(i2)=x1.x2(i1);
    else
        x2=mergeStruct(x2,trimStruct(x1,i1));
        i2=i2+1;
    end
end

    
%%
function test1
clear
clf
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/CN/CCF')
X=load_tsv('REBC_primary_pair_QC_25Jan2020_383.alleliccapseg_pp_ccf_fit_v3.tsv')
X.pair_id=X.sample;
pid1=unique(sort(X.sample));
%tab(A.sample)
tab(X.pair_id)

pid=sort(unique(X.pair_id));
N=length(pid);
db=0.0025;
acrb=0:db:5;
acrr=-0.5:db:0.5;

aL=num2chrom(1:23);
xL=xhg19(aL,zeros(size(aL)));
xM=round(conv(xL,[1 1]/2, 'same')); xM(end)=[];
    
X.x1=xhg19(X.Chromosome,X.Start_bp);
X.x2=xhg19(X.Chromosome,X.End_bp);

X2=[]
for i=1:N
    x1=trimStruct(X,ismember(X.pair_id,pid(i)));
    [~,k]=sort(x1.x1);
    x1=trimStruct(x1,k);
    x2=collapse_AllelicCapSeg_PP_CCF_fit_v3(x1);
    if isempty(X2)
        X2=x2;
    else
        X2=mergeStruct(X2,x2);
    end
    [i x1.N x2.N X2.N]
end

t1=tab(X.pair_id)
t2=tab(X2.pair_id)
[~,k]=sort(t1.n)
t1=trimStruct(t1,k)
[i m]=ismember(t1.x,t2.x)
semilogy([t1.n t2.n(m)])

%%
h1=subplot(1,1,1)
h1.Box='on'
for i=34:N
    clf
    pid(i)
    x1=trimStruct(X,ismember(X.pair_id,pid(i)));
    x2=trimStruct(X2,ismember(X2.pair_id,pid(i)));
    if sum(x2.IS_SCNA.*x2.length)<10e6,continue; end
    pid(i)
    if ~sum((x2.IS_SCNA.*(x2.length>10e6).*ismember(x2.Chromosome,22))),continue; end
    pid(i)
    k=find(ismember(x2.Chromosome,22));
    x22=trimStruct(x2,k);
    CCF=max(x22.NA_CCF);        
    line([xL xL]',[0+0*xL 4+0*xL]','color',0.75*[1 1 1],'linestyle','--');
    set(gca,'Box','on')
    for i1=1:x1.N
        xx=[x1.x1(i1) x1.x1(i1)+x1.length(i1) x1.x1(i1)+x1.length(i1) x1.x1(i1) x1.x1(i1)];
        yy=[x1.mu_minor(i1)-x1.sigma_minor(i1)/1  x1.mu_minor(i1)-x1.sigma_minor(i1)/1 x1.mu_minor(i1)+x1.sigma_minor(i1)/1 x1.mu_minor(i1)+x1.sigma_minor(i1)/1 x1.mu_minor(i1)-x1.sigma_minor(i1)/1];
        patch(xx,yy,'blue','EdgeColor','none','FaceAlpha',0.5);
        yy=[x1.mu_major(i1)-x1.sigma_major(i1)/1  x1.mu_major(i1)-x1.sigma_major(i1)/1 x1.mu_major(i1)+x1.sigma_major(i1)/1 x1.mu_major(i1)+x1.sigma_major(i1)/1 x1.mu_major(i1)-x1.sigma_major(i1)/1];
        patch(xx,yy,'red','EdgeColor','none','FaceAlpha',0.5);
        line([x1.x1(i1) x1.x2(i1)],x1.cn_minor_cluster(i1)*[ 1 1]+1e-3,'linestyle','-','Color','b');
        line([x1.x1(i1) x1.x2(i1)],x1.cn_major_cluster(i1)*[ 1 1],'linestyle','-','Color','r');
        set(gca,'xtick',xM,'xticklabel',aL,'ylim',[0 4]);
    end
    for i2=1:x2.N
        xx=[x2.x1(i2) x2.x1(i2)+x2.length(i2) x2.x1(i2)+x2.length(i2) x2.x1(i2) x2.x1(i2)];
        yy=[x2.mu_minor(i2)-x2.sigma_minor(i2)/1  x2.mu_minor(i2)-x2.sigma_minor(i2)/1 x2.mu_minor(i2)+x2.sigma_minor(i2)/1 x2.mu_minor(i2)+x2.sigma_minor(i2)/1 x2.mu_minor(i2)-x2.sigma_minor(i2)/1];
        patch(xx,yy+0.1,'b','EdgeColor','none','FaceAlpha',0.5);
        yy=[x2.mu_major(i2)-x2.sigma_major(i2)/1  x2.mu_major(i2)-x2.sigma_major(i2)/1 x2.mu_major(i2)+x2.sigma_major(i2)/1 x2.mu_major(i2)+x2.sigma_major(i2)/1 x2.mu_major(i2)-x2.sigma_major(i2)/1];
        patch(xx,yy+0.1,'r','EdgeColor','none','FaceAlpha',0.5);
        line([x2.x1(i2) x2.x2(i2)],x2.cn_minor_cluster(i2)*[ 1 1]+1e-3,'linestyle','-','Color','b');
        line([x2.x1(i2) x2.x2(i2)],x2.cn_major_cluster(i2)*[ 1 1],'linestyle','-','Color','r');
        set(gca,'xtick',xM,'xticklabel',aL,'ylim',[0 4]);
    end
    title(sprintf('%s pur=%.2f ploidy=%.2f N1=%d N2=%d CCF=%.2f',pid{i},x1.purity(1),x1.ploidy(1),x1.N,x2.N,CCF))
    ylim([0 2.5])
    xlim([0 3e9])
    xlabel('chromosome')
    ylabel('allelic copy ratio')
    pause()
end

%%
%22q DEL clonal 


