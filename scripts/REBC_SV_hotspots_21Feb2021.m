clear
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/')
%PAIR_FCDATE='11Apr2020';
%SAMPLE_FCDATE='25Sep2019';

% Three REBC pair_sets - not including mets or second primaries
PAIR_SETS0={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};

% filter original REBC pair sets for hotspots
AREA1='SV/'
PAIR_SETS1=PAIR_SETS0;

% filter met and second primary REBC pair sets for hotspots
AREA1='mets/SV/'
PAIR_SETS1={'REBC_met_concurrent_primary_pair_60','REBC_met_concurrent_corresponding_pair_59'};


X=[]
for PAIR_SET=PAIR_SETS0
    PAIR_SET
    % load SV lists attribute SV_BP_CCF_v3 from task tsvcat_SV_BP_CCF_v3  
    X1=load_tsv(['SV/' PAIR_SET{1} '.SV_BP_CCF_v3.tsv']);
    X1.PAIR_SET=repmat(PAIR_SET,size(X1.chr1));
    if isempty(X)
        X=X1;
    else
        X=mergeStruct(X,X1)
    end
end

str1=1-X.str1*2;
str2=1-X.str2*2;
X.strand1=repmat({'+'},size(X.str1))
k=find(X.str1>0);
X.strand1(k)={'-'};
X.strand2=repmat({'+'},size(X.str1))
k=find(X.str2>0);
X.strand2(k)={'-'};

% generate HOTSPOT list of breaks from these three pairsets

X.eid=regexprep(cellstr((strcat(num2str(X.chr1),':',num2str(X.pos1),':',X.strand1,',',num2str(X.chr2),':',num2str(X.pos2),':',X.strand2,',',X.gene1,'-',X.gene2))),' ','')
tab(X.eid)
t=tab(X.eid)
HOTSPOTS=t.x(t.n>1)

%fragile sites
% PCAWG Nature. 2020; 578(7793): 112?121.
% Published online 2020 Feb 5. doi: 10.1038/s41586-019-1913-9
 
F=load_tsv('SV/FragileSites/PCAWG_FragileSites_STable5_LiEtAl_Nature_2020.tsv')
F.x1=xhg19(F.chrom,F.start);
F.x2=xhg19(F.chrom,F.xEnd);


for PAIR_SET=PAIR_SETS1
    PAIR_SET
    % load SV lists attribute SV_BP_CCF_v3 from task tsvcat_SV_BP_CCF_v3  
    X=load_tsv([AREA1 PAIR_SET{1} '.SV_BP_CCF_v3.tsv']);
    X.PAIR_SET=repmat(PAIR_SET,size(X.chr1));   

    str1=1-X.str1*2;
    str2=1-X.str2*2;
    X.strand1=repmat({'+'},size(X.str1))
    k=find(X.str1>0);
    X.strand1(k)={'-'};
    X.strand2=repmat({'+'},size(X.str1))
    k=find(X.str2>0);
    X.strand2(k)={'-'};
    
    X.eid=regexprep(cellstr((strcat(num2str(X.chr1),':',num2str(X.pos1),':',X.strand1,',',num2str(X.chr2),':',num2str(X.pos2),':',X.strand2,',',X.gene1,'-',X.gene2))),' ','')


    X.fragileSitePCAWG=repmat({''},size(X.chr1))
    [i1,m1] = point_overlap_intervals(X.x1,F);
    tab(i1)
    X.fragileSitePCAWG(i1)=F.CFS(m1(i1));
    [i1,m1] = point_overlap_intervals(X.x2,F);
    tab(i1)
    X.fragileSitePCAWG(i1)=F.CFS(m1(i1));

    X1=trimStruct(X,~ismember(X.eid,HOTSPOTS))
    
    % SV list
    printStruct(X1,-1,[AREA1 PAIR_SET{1} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv']);
end

% load as attribute back to Firecloud 
mkdir('/tmp/SV')
fid=fopen('/tmp/SV/load.tsv','w')
fprintf(fid,'update_pair_id	SV_BP_CCF_v3_nohot_tsv\n')

for PAIR_SET=PAIR_SETS1
    PAIR_SET
    PS=load_tsv(['FC/' PAIR_SET{1} '.pair_set.tsv']);
    %P=load_tsv(['FC/REBC.pairs.' PAIR_FCDATE '.tsv']);
    X=load_tsv([AREA1 PAIR_SET{1} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv'])
    N=length(PS.pair_id)
    for i=1:N
        x1=trimStruct(X,ismember(X.individual,PS.pair_id(i)));
        [PS.pair_id{i} '   ' num2str(x1.N)]
        file0=['/tmp/SV/' PS.pair_id{i} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv'];
        file1=['gs://fc-035f5652-acf7-4642-abb7-e8c10848c8ed/' AREA1 PS.pair_id{i} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv'];
        printStruct(x1,-1,file0);
        cmd=['gsutil cp ' file0 '  ' file1];
        unix(cmd)
        fprintf(fid,'%s %s\n',PS.pair_id{i},file1);        
    end
end
fclose(fid)
cmd('bbedit /tmp/SV/load.tsv')