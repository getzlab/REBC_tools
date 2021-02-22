clear
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};

AREA='mets/SV/'
PAIR_SETS={'REBC_met_concurrent_primary_pair_60','REBC_met_concurrent_corresponding_pair_59'};
NPS=length(PAIR_SETS)
for i=1:NPS
    PSET=PAIR_SETS{i}
    unix(['ls -latr SV/' PSET '*.tsv'])
    
    %X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filtNA01.24May2020.tsv'])
    X=load_tsv([AREA PSET '.SV_BP_CCF_v3.noHot.21Feb2021.tsv'])
    %X0=load_tsv([AREA PSET '.aggregated.all.tsv'])
    X0=load_tsv([AREA PSET '.SV_breakpointer_fix_sample.tsv'])
    
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
    
    printStruct(X2,-1,[AREA PSET '.aggregated.all.clustered.' TODAY '.tsv'])
    %XX=trimStruct(X2,X2.passcluster&(~X2.pass))
    
    
    X3=trimStruct(X2,X2.pass|X2.passcluster)
    printStruct(X2,find(X2.pass|X2.passcluster),[AREA PSET '.SV_BP_CCF_v3.noHot.passcluster.' TODAY '.tsv'])
    
    % printStruct(X3,ismember(X3.individual,'REBC-ACCF-TP-NT')&ismember(X3.gene1,'PAX8'))
    % printStruct(X3,ismember(X3.individual,'REBC-AC9K-TP-NT')&ismember(X3.gene1,'RET'))
end

%%

% load as attribute back to Firecloud 
mkdir('/tmp/SV')
fid=fopen('/tmp/SV/load.tsv','w')
fprintf(fid,'update_pair_id	SV_BP_CCF_v3_nohot_passcluster_tsv\n')
DAY='21Feb2021'
for PAIR_SET=PAIR_SETS
    PAIR_SET
    PS=load_tsv(['FC/' PAIR_SET{1} '.pair_set.tsv']);
    %P=load_tsv(['FC/REBC.pairs.' PAIR_FCDATE '.tsv']);
    X=load_tsv([AREA PAIR_SET{1} '.SV_BP_CCF_v3.noHot.passcluster.' DAY '.tsv'])
    N=length(PS.pair_id)
    for i=1:N
        x1=trimStruct(X,ismember(X.individual,PS.pair_id(i)));
        [PS.pair_id{i} '   ' num2str(x1.N)]
        file0=['/tmp/SV/' PS.pair_id{i} '.SV_BP_CCF_v3.noHot.passcluster.' DAY '.tsv'];
        file1=['gs://fc-035f5652-acf7-4642-abb7-e8c10848c8ed/' AREA PS.pair_id{i} '.SV_BP_CCF_v3.noHot.passcluster.' DAY '.tsv'];
        printStruct(x1,-1,file0);
        cmd=['gsutil cp ' file0 '  ' file1];
        unix(cmd)
        fprintf(fid,'%s %s\n',PS.pair_id{i},file1);        
    end
end
fclose(fid)
cmd('bbedit /tmp/SV/load.tsv')

%%clear
X3=load_tsv('/Volumes/GoogleDrive/My Drive/Cancer/REBC/mets/SV/REBC_met_concurrent_primary_pair_60.SV_BP_CCF_v3.noHot.passcluster.21Feb2021.tsv')
 
FUSION_DRIVERS={ 'RET', 'NTRK3',  'NTRK1', 'BRAF','ALK', 'THADA', 'LTK', 'IGF2', 'IGF2BP3','PPARG'};
XD=trimStruct(X3,ismember(X3.gene1,FUSION_DRIVERS)|ismember(X3.gene2,FUSION_DRIVERS))


tab(XD.ID)
printStruct(XD,ismember(XD.individual,'REBC-ACC9-TP-NB'))

