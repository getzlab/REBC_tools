clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
for i=1:3
    PSET=PAIR_SETS{i}
    unix(['ls -latr SV/' PSET '*.tsv'])
    
    %X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filt.24May2020.tsv'])
    %X0=load_tsv(['SV/' PSET '.aggregated.all.clustered.26May2020.tsv'])
    X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.tsv'])
    %X0=load_tsv(['SV/' PSET '.aggregated.all.filtNA01.24May2020.tsv'])
    X0=load_tsv(['SV/' PSET '.aggregated.all.tsv'])
    
    X0.NALG=(X0.dRanger>0)+(X0.pcawg_snowman>0)+(X0.SvABA>0)+(X0.Manta>0);
    if isfield(X,'balanced')
        X=RenameField(X,'balanced','balanced0')
    end
    %
    % remove SvABA pcawg_snowman_only
    X0=trimStruct(X0,~((X0.NALG==1)&(X0.SvABA|X0.pcawg_snowman)))
    %
    STR='+-'
    str1=cellstr(STR(X.str1+1)');
    str2=cellstr(STR(X.str2+1)');
    X.eid=regexprep(strcat(X.individual,':',num2str(X.chr1),':',num2str(X.pos1),':',str1,',',num2str(X.chr2),':',num2str(X.pos2),':',str2,',',X.gene1,'-',X.gene2),' ','');
    str1=cellstr(STR(X0.str1+1)');
    str2=cellstr(STR(X0.str2+1)');
    X0.eid=regexprep(strcat(X0.individual,':',num2str(X0.chr1),':',num2str(X0.pos1),':',str1,',',num2str(X0.chr2),':',num2str(X0.pos2),':',str2,',',X0.gene1,'-',X0.gene2),' ','');
    
    % SV calls that failed filtering for SV_BP_CCF_v3.noHot.filtNA01.24May2020
    [i1 m1]=ismember(X0.eid,X.eid);
    XX=trimStruct(X0,find(~i1));
    
    X1=mergeStruct(X,XX);
    X1.pass=ismember(X1.eid,X.eid);
    [~,k]=sort(X1.eid);
    X1=trimStruct(X1,k);
    X1.x1=xhg19(X1.chr1,X1.pos1);
    X1.x2=xhg19(X1.chr2,X1.pos2);
    %
    WB=1000
    % window beefed up to 50kb
    WC=50000
    [X2,Z2] = SV_cluster(X1,WB,WC)
%     
%     if(0)
%   %      printStruct(X,ismember(X.individual,'REBC-ACCF-TP-NT')&ismember(X.gene1,'PAX8'))
%   %      printStruct(X1,ismember(X1.individual,'REBC-ACCF-TP-NT')&ismember(X1.gene1,'PAX8'))
%     end
%     
%     
    % match clustered SVs to filtered SVs
    %X1.pass=ismember(X1.eid,X.eid);
    % SV is pass or within a balanced,inv2, or complex
    %trimStruct(X1,X1.pass|(X1.bal>0)|(X1.inv2>0)|(X1.complex>0))
    tab(X2.individual)
    id=unique(X2.individual)
    % SV is pass of part of a balanced,inv2, or complex with a pass
    X2.passcluster=0*X2.chr1;
    for id1=id'
        id1
        k=find(ismember(X2.individual,id1));
        ub=unique(X2.balanced(k)); ub=ub(ub>0);
        for ub1=ub'
            kb=k(find(ismember(X2.balanced(k),ub1)));
            X2.passcluster(kb)=any(X2.pass(kb));
        end
        ui=unique(X2.inversion2(k)); ui=ui(ui>0);
        for ui1=ui'
            ki=k(find(ismember(X2.inversion2(k),ui1)));
            X2.passcluster(ki)=any(X2.pass(ki));
        end
        uc=unique(X2.complex(k)); uc=uc(uc>0);
        for uc1=uc'
            kc=k(find(ismember(X2.complex(k),uc1)));
            X2.passcluster(kc)=any(X2.pass(kc));
        end
    end
    
    printStruct(X2,-1,['SV/' PSET '.aggregated.all.clustered.' TODAY '.tsv'])
    %XX=trimStruct(X2,X2.passcluster&(~X2.pass))
    
    
    X3=trimStruct(X2,X2.pass|X2.passcluster)
    printStruct(X2,find(X2.pass|X2.passcluster),['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.passcluster.' TODAY '.tsv'])
    
    % printStruct(X3,ismember(X3.individual,'REBC-ACCF-TP-NT')&ismember(X3.gene1,'PAX8'))
    % printStruct(X3,ismember(X3.individual,'REBC-AC9K-TP-NT')&ismember(X3.gene1,'RET'))
end

%%
clear
X3=load_tsv('/GoogleDrive/Cancer/REBC/SV/REBC_primary_pair_393.SV_BP_CCF_v3.noHot.filtNA01.24May2020.passcluster.28Jun2020.tsv')
printStruct(X3,ismember(X3.individual,'REBC-ACCF-TP-NT')&ismember(X3.gene1,'PAX8'))
 
FUSION_DRIVERS={ 'RET', 'NTRK3',  'NTRK1', 'BRAF','ALK', 'THADA', 'LTK', 'IGF2', 'IGF2BP3','PPARG'};
XD=trimStruct(X3,ismember(X3.gene1,FUSION_DRIVERS)|ismember(X3.gene2,FUSION_DRIVERS))

tab(XD.ID)
printStruct(XD,ismember(XD.individual,'REBC-ACC9-TP-NB'))
% 

%%
clear
cd('/GoogleDrive/Cancer/REBC/')
FUSION_DRIVERS={ 'RET', 'NTRK3',  'NTRK1', 'BRAF','ALK', 'THADA', 'LTK', 'IGF2', 'IGF2BP3','PPARG'};
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
MTDAY='16Jun2020'; VERSION='V12'
M0=load_tsv(['MainTable/REBC_THCA_TABLE_' VERSION '.' MTDAY '.tsv']);
M0=trimStruct(M0,~ismember(M0.WGS_TP,'n/a'))
M0.N
SVDAY='28Jun2020'
for i=1:1
    PSET=PAIR_SETS{i}
    if PSET(1)=='R'
        M=trimStruct(M0,cellfun(@length,regexp(M0.LABELED_REBC_ID,'^R'))>0)
    else
        M=trimStruct(M0,cellfun(@length,regexp(M0.LABELED_REBC_ID,'^T'))>0)
    end
    M=rmfield(M,'N');
    
    X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.passcluster.' SVDAY '.tsv'])
    X.LABELED_REBC_ID=cellfun(@(x) x(1), regexp(X.individual,'-TP','split'))
    X=trimStruct(X,ismember(X.LABELED_REBC_ID,M.LABELED_REBC_ID))
    ID=M.LABELED_REBC_ID;
    Z=[]
    Z.ID=ID;
    [i1 m1]=ismember(ID,M.LABELED_REBC_ID)
    Z.DOSE0712_STO_ARITH_MEAN_THY(i1,1)=M.DOSE0712_STO_ARITH_MEAN_THY(m1(i1))
    Z.IREP2_POC_50TH(i1,1)=M.IREP2_POC_50TH(m1(i1))
    Z.RE2_AGE_SURG(i1,1)=M.RE2_AGE_SURG(m1(i1))
    Z.RE2_AGE_EXP(i1,1)=M.RE2_AGE_EXP(m1(i1))
    Z.WGS_TP_N_QCPASS(i1,1)=M.WGS_TP_N_QCPASS(m1(i1))
    
    % DOSE bins
    DB=[0,1,100,200,500,1000,1e8]'
    DBL={'     0','    1-99','   100-199','  200-499', ' 500-1000','1000+'}'
    xd=str2double(Z.DOSE0712_STO_ARITH_MEAN_THY); nanmin(xd); xd(isnan(xd))=0
    [n,bd]=histc(xd,DB-0.00001)
    Z.DOSE0712_STO_ARITH_MEAN_THY_BIN=DBL(bd)
    Z.DOSE0712_STO_ARITH_MEAN_THY=str2double(Z.DOSE0712_STO_ARITH_MEAN_THY)
    Z.DOSE0712_STO_ARITH_MEAN_THY(isnan(Z.DOSE0712_STO_ARITH_MEAN_THY))=0
    
    % check POC bins
    PB=[0,1,25,50,75,1e8]
    PBL={'    0','   1-24','  25-49',' 50-74','75+'}'
    xp=str2double(Z.IREP2_POC_50TH); nanmin(xp)
    xp(isnan(xp))=0
    Z.IREP2_POC_50TH=str2double(Z.IREP2_POC_50TH)
    Z.IREP2_POC_50TH(isnan(Z.IREP2_POC_50TH))=0
    [n1,bp]=histc(xp,PB-0.00001)
    Z.IREP2_POC_50TH_BIN=PBL(bp)
    
    
    Z.RE2_LATENCY=Z.RE2_AGE_SURG- Z.RE2_AGE_EXP;
    for i=1:length(ID)
        x1=trimStruct(X,ismember(X.LABELED_REBC_ID,Z.ID(i)));
        Z.nSV(i,1)=x1.N;
        Z.nSVbalanced(i,1)=sum(x1.balanced>0);
        Z.nSVinversion2(i,1)=sum(x1.inversion2>0);
        Z.nSVcomplex(i,1)=sum(x1.complex>0);
        Z.nSVbalancedNOTcomplex(i,1)=sum((x1.balanced>0)&(x1.complex==0));
        Z.nSVinv2NOTcomplex(i,1)=sum((x1.inversion2>0)&(x1.complex==0));
        Z.nSVisolated(i,1)=sum((x1.balanced+x1.inversion2+x1.complex)==0);
        % mean distance for balanced 
        
        % cluster counts
        kc=find(x1.complex>0);
        Z.nComplex(i,1)=length(unique(x1.complex(kc)));
        Z.distComplex(i,1)=median([x1.clusterDeltaPos1(kc); x1.clusterDeltaPos2(kc)]); 
        kb=find((x1.complex==0)&(x1.balanced>0));
        Z.nBalanced(i,1)=length(unique(x1.balanced(kb)));
        Z.distBalanced(i,1)=median(abs([x1.clusterDeltaPos1(kb); x1.clusterDeltaPos2(kb)])); 
        ki=find((x1.complex==0)&(x1.inversion2>0));
        Z.nInversion2(i,1)=length(unique(x1.inversion2(ki)));
        Z.distInversion2(i,1)=median(abs([x1.clusterDeltaPos1(ki); x1.clusterDeltaPos2(ki)])); 
        
        % complex occupancy
        kc=unique(x1.complex); kc=kc(kc>0);
        n1=[];
        for c1=kc(:)'
            c1;
            n1=[n1 sum(x1.complex==c1)];
        end
        Z.complexMeanSV(i,1)=0;
        Z.complexMedianSV(i,1)=0;
        Z.complexMaxSV(i,1)=0;
        if length(n1)>0
            n1;
            Z.complexMeanSV(i,1)=mean(n1);
            Z.complexMedianSV(i,1)=median(n1);
            Z.complexMaxSV(i,1)=max(n1);
        end
        
        % drivers 
        xd=trimStruct(x1,ismember(x1.gene1,FUSION_DRIVERS)|ismember(x1.gene2,FUSION_DRIVERS));
        Z.nSVdriver(i,1)=xd.N;
        Z.nSVbalancedDriver(i,1)=sum(xd.balanced>0);
        Z.nSVInversion2Driver(i,1)=sum(xd.inversion2>0);
        Z.nSVcomplexDriver(i,1)=sum(xd.complex>0);
        Z.nSVbalancedNOTcomplexDriver(i,1)=sum((xd.balanced>0)&(xd.complex==0));
        Z.nSVinv2NOTcomplexDriver(i,1)=sum((xd.inversion2>0)&(xd.complex==0));
        Z.nSVisolatedDriver(i,1)=sum((xd.balanced+xd.inversion2+xd.complex)==0);
        %Z.nSVsimpleDriver(i,1)=Z.nSVisolatedDriver(i,1)+Z.nSVbalancedNOTcomplexDriver(i,1)+Z.nSVinv2NOTcomplexDriver(i,1);
        
        % driver cluster counts
        kc=find(xd.complex>0);
        Z.nComplexDriver(i,1)=length(unique(xd.complex(kc)));
        Z.distComplexDriver(i,1)=median([xd.clusterDeltaPos1(kc); xd.clusterDeltaPos2(kc)]);      
        kb=find((xd.complex==0)&(xd.balanced>0));
        Z.nBalancedDriver(i,1)=length(unique(xd.balanced(kb)));
        Z.distBalancedDriver(i,1)=median(abs([xd.clusterDeltaPos1(kb); xd.clusterDeltaPos2(kb)]));              
        ki=find((xd.complex==0)&(xd.inversion2>0));
        Z.nInversion2Driver(i,1)=length(unique(xd.inversion2(ki)));
        Z.distInversion2Driver(i,1)=median(abs([xd.clusterDeltaPos1(ki); xd.clusterDeltaPos2(ki)])); 
     
        % complex occupancy
        kc=unique(xd.complex); kc=kc(kc>0);
        n1=[];
        for c1=kc(:)'
            c1;
            n1=[n1 sum(xd.complex==c1)];
        end
        Z.complexMeanSVDriver(i,1)=0;
        Z.complexMedianSVDriver(i,1)=0;
        Z.complexMaxSVDriver(i,1)=0;
        if length(n1)>0
            n1;
            Z.complexMeanSVDriver(i,1)=mean(n1);
            Z.complexMedianSVDriver(i,1)=median(n1);
            Z.complexMaxSVDriver(i,1)=max(n1);
        end

    
    end
    
    Z.nTot=Z.nSVisolated+Z.nComplex+ Z.nBalanced+ Z.nInversion2;
    
    printStruct(Z)
    printStruct(Z,-1,['SV/' PSET '.passcluster.' SVDAY '.clustering_metrics.' TODAY '.tsv'])
    
end

%%
clear
Z=load_tsv('SV/REBC_primary_pair_393.passcluster.28Jun2020.clustering_metrics.30Jun2020.tsv')
 

% driver that looks like radiation (not complex - balanced or inverse2 with dist<20?

myhist(Z.distBalancedDriver(Z.nBalancedDriver>0),0:50)
myhist(Z.distInversion2Driver(Z.nInversion2Driver>0),0:50)

db=unique(Z.DOSE0712_STO_ARITH_MEAN_THY_BIN)
D=[]
D.db=db
ndb=length(db)

histstr(Z.DOSE0712_STO_ARITH_MEAN_THY_BIN,db)

Z1=trimStruct(Z,Z.WGS_TP_N_QCPASS>0)
Z1.dist=nanmin([Z1.distBalancedDriver Z1.distInversion2Driver],[],2)

for i=1:ndb
    z1=trimStruct(Z1,ismember(Z1.DOSE0712_STO_ARITH_MEAN_THY_BIN,D.db(i)))
    z1=trimStruct(z1,z1.nSVdriver>0)
    D.DOSE0712_STO_ARITH_MEAN_THY(i,1)=mean((z1.DOSE0712_STO_ARITH_MEAN_THY));
    D.SVisolatedDriver(i,1)=sum((z1.nSVisolatedDriver>0)&(z1.nComplexDriver==0)&(z1.nBalancedDriver==0)&(z1.nInversion2Driver==0));
    D.SVbalancedDriver(i,1)=sum(z1.nBalancedDriver>0);
    D.SVInversion2Driver(i,1)=sum(z1.nInversion2Driver>0);
    D.SVComplexDriver(i,1)=sum(z1.nComplexDriver>0);
    
    D.unknown(i,1)=sum((z1.nSVisolatedDriver>0)&(z1.nComplexDriver==0)&(z1.nBalancedDriver==0)&(z1.nInversion2Driver==0));
    D.NHEJL(i,1)=sum((z1.nComplexDriver==0)&((z1.nBalancedDriver>0)|(z1.nInversion2Driver>0))&(z1.dist<15));
    D.NHEJH(i,1)=sum((z1.nComplexDriver==0)&((z1.nBalancedDriver>0)|(z1.nInversion2Driver>0))&(z1.dist>=15));
    
end


[~,k]=sort(D.DOSE0712_STO_ARITH_MEAN_THY)
D=trimStruct(D,k)

subplot(2,1,1)
bar([D.SVisolatedDriver D.SVbalancedDriver D.SVInversion2Driver D.SVComplexDriver],'stacked')
lab={'isolatedDriver','balancedDriver','Inversion2Driver','ComplexDriver'}
legend(lab)
set(gca,'xticklabel',D.db)
xlim([0.5 6.5])
grid on
xlabel('DOSE bin (mGy)'); ylabel('tumors')
subplot(2,1,2)
bar([D.SVisolatedDriver D.SVbalancedDriver D.SVInversion2Driver D.SVComplexDriver]./repmat(sum([D.SVisolatedDriver D.SVbalancedDriver D.SVInversion2Driver D.SVComplexDriver],2),[1,4]),'stacked')
lab={'isolatedDriver','balancedDriver','Inversion2Driver','ComplexDriver'}
legend(lab)
set(gca,'xticklabel',D.db)
xlim([0.5 6.5])
grid on
xlabel('DOSE bin (mGy)'); ylabel('relative proportions of tumors')
  

figure(2)
subplot(2,1,1)
h=bar([D.unknown D.NHEJL D.NHEJH D.SVComplexDriver],'stacked')
lab={'unknown','NHEJ d<15bp Driver','NHEJ d>=15bp Driver','Complex Driver'}
h(1).FaceColor=0.5*[1 1 1]
legend(lab)
set(gca,'xticklabel',D.db)
xlim([0.5 6.5])
grid on
xlabel('DOSE bin (mGy)'); ylabel('tumors')
subplot(2,1,2)
h=bar([D.unknown D.NHEJL D.NHEJH D.SVComplexDriver]./repmat(sum([D.unknown D.NHEJL D.NHEJH D.SVComplexDriver],2),[1,4]),'stacked')
lab={'unknown','NHEJ d<15bp Driver','NHEJ d>=15bp Driver','Complex Driver'}
h(1).FaceColor=0.5*[1 1 1]
legend(lab,'location','E')
set(gca,'xticklabel',D.db)
xlim([0.5 6.5])
grid on
xlabel('DOSE bin (mGy)'); ylabel('relative proportions of tumors')
  
%%
clear
cd('/GoogleDrive/Cancer/REBC/')
FUSION_DRIVERS={ 'RET', 'NTRK3',  'NTRK1', 'BRAF','ALK', 'THADA', 'LTK', 'IGF2', 'IGF2BP3','PPARG'};
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
MTDAY='16Jun2020'; VERSION='V12'
M0=load_tsv(['MainTable/REBC_THCA_TABLE_' VERSION '.' MTDAY '.tsv']);
M0=trimStruct(M0,~ismember(M0.WGS_TP,'n/a'))
M0.N
SVDAY='28Jun2020'
PSET=PAIR_SETS{1}
X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.passcluster.' SVDAY '.tsv'])
  
M=trimStruct(M0,cellfun(@length,regexp(M0.LABELED_REBC_ID,'^R'))>0)
M=trimStruct(M,M.WGS_TP_N_QCPASS>0)
X.LABELED_REBC_ID=cellfun(@(x) x(1), regexp(X.individual,'-TP','split'))
X=trimStruct(X,ismember(X.LABELED_REBC_ID,M.LABELED_REBC_ID))
  
    

figure(3)
clf
subplot(2,1,1)
x=[X.clusterDeltaPos1((X.balanced>0)|(X.inversion2>0)); X.clusterDeltaPos2((X.balanced>0)|(X.inversion2>0))];
myhist(x,-200:2:200)
xlabel('distance between reciprocal breaks (bp)')
ylabel('balanced or inversion2 SV breaks ')
FUSION_DRIVERS={ 'RET', 'NTRK3',  'NTRK1', 'BRAF','ALK', 'THADA', 'LTK', 'IGF2', 'IGF2BP3','PPARG'};
XD=trimStruct(X,ismember(X.gene1,FUSION_DRIVERS)|ismember(X.gene2,FUSION_DRIVERS));
subplot(2,1,2)
xd=[XD.clusterDeltaPos1((XD.balanced>0)|(XD.inversion2>0)); XD.clusterDeltaPos2((XD.balanced>0)|(XD.inversion2>0))];
myhist(xd,-200:2:200)
xlabel('distance between reciprocal breaks (bp)')
ylabel({'Driver',' balanced or inversion2 SV breaks '})


% subplot(2,1,1)
% myhist(X.clusterDeltaPos1((X.balanced>0)|(X.inversion2>0)),-200:2:200)
% hold on;
% myhist(X.clusterDeltaPos2((X.balanced>0)|(X.inversion2>0)),-200:2:200)



    %%
clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
for i=1:3
    PSET=PAIR_SETS{i}
    unix(['ls -latr SV/' PSET '*.tsv'])

    X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.tsv'])

    WB=1000
    WC=50000
    [X1,Z] = SV_cluster(X,WB,WC)
    printStruct(X1,-1,['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.cluster.50kb.' TODAY '.tsv'])
    printStruct(Z,-1,['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.cluster.50kb.summary.' TODAY '.tsv'])
end

%%
clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
PSET=PAIR_SETS{1}
unix(['ls -latr SV/' PSET '*.tsv'])
DAY= '14Jun2020'
X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.04May2020.pcawg_clustered.26May2020.cluster.50kb' DAY '.tsv'])

figure(1)
clf
P=[]
P.log=1
ALGS={'dRanger','SvABA','pcawg_snowman','Manta'}
h=upset_plot(X,ALGS,P)
title(h(1),PSET)
f=['SV/plots/' PSET '.SV_raw_UpSet.SV_BP_CCF_v3.noHot.04May2020.pcawg_clustered.26May2020.cluster.50kb.' TODAY ]
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')
saveas(gcf,[f '.jpg'],'jpeg')

k=find((X.dRanger<1)&(X.Manta<1))
printStruct(X,k)

%%
clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
PSET=PAIR_SETS{1}
unix(['ls -latr SV/' PSET '*.tsv'])
DAY= '14Jun2020'
fSV=['SV/' PSET '.aggregated.all.clustered.passpart.' DAY '.tsv']
X=load_tsv(fSV)

figure(1)
clf
P=[]
P.log=1
ALGS={'dRanger','SvABA','pcawg_snowman','Manta'}
h=upset_plot(X,ALGS,P)
t=regexprep(fSV,'SV/',''); t=regexprep(t,'.tsv','')
title(h(1),t)
f=[ 'SV/plots/' t ]
text(h(1),0.7,0.9,sprintf('total #SV: %d',length(X.chr1)),'units','normalized','fontsize',14) 
saveas(gcf,[f '.svg'],'svg')
saveas(gcf,[f '.eps'],'epsc')
saveas(gcf,[f '.jpg'],'jpeg')

k=find((X.dRanger<1)&(X.Manta<1))
printStruct(X,k)


printStruct(X,X.passpart&~X.pass)

%% 
tab(X.ID)


