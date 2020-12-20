clear
PAIR_FCDATE='11Apr2020';
SAMPLE_FCDATE='25Sep2019';
MAFDATE='07Apr2020';
SV_CLONAL_CCF_THRESHOLD=0.75;
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
cd('/Volumes/GoogleDrive/My Drive/Cancer/REBC/')

X=[]
for PAIR_SET=PAIR_SETS
    PAIR_SET
    % SV 
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


X.eid=regexprep(cellstr((strcat(num2str(X.chr1),':',num2str(X.pos1),':',X.strand1,',',num2str(X.chr2),':',num2str(X.pos2),':',X.strand2,',',X.gene1,'-',X.gene2))),' ','')
tab(X.eid)
t=tab(X.eid)
hotspot=t.x(t.n>1)
[i m]=ismember(X.eid,t.x);
tab(i)
X.recurrence(i,1)=t.n(m(i))-1;


k=find(ismember(X.PAIR_SET,PAIR_SETS(1)))
tab(X.eid(k))
t1=tab(X.eid(k))
hotspot1=t1.x(t1.n>1)

hotspotx=hotspot(find(~ismember(hotspot,hotspot1)))
XX=trimStruct(X,ismember(X.eid,hotspotx))

% manual review all XX were artifacts with one germline 
% 2:43374382:-,2:65045476:+,ZFP36L2-SLC1A4 in REBC-AC92

tab(XX.individual)

XH=trimStruct(X,ismember(X.eid,hotspot))
[~,k]=sort(XH.eid)
XH=printStruct(XH,k)

%fragile sites
% PCAWG Nature. 2020; 578(7793): 112?121.
% Published online 2020 Feb 5. doi: 10.1038/s41586-019-1913-9
 
F=load_tsv('SV/FragileSites/PCAWG_FragileSites_STable5_LiEtAl_Nature_2020.tsv')
F.x1=xhg19(F.chrom,F.start);
F.x2=xhg19(F.chrom,F.xEnd);

X.fragileSitePCAWG=repmat({''},size(X.chr1))
[i1,m1] = point_overlap_intervals(X.x1,F);
tab(i1)
X.fragileSitePCAWG(i1)=F.CFS(m1(i1));
[i1,m1] = point_overlap_intervals(X.x2,F);
tab(i1)
X.fragileSitePCAWG(i1)=F.CFS(m1(i1));

XPASS=trimStruct(X,~ismember(X.eid,hotspot))


for PAIR_SET=PAIR_SETS
    PAIR_SET
    % SV 
    X1=trimStruct(XPASS,ismember(XPASS.PAIR_SET,PAIR_SET))
    printStruct(X1,-1,['SV/' PAIR_SET{1} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv']);
end

mkdir('/tmp/SV')
fid=fopen('/tmp/SV/load.tsv','w')
fprintf(fid,'update_pair_id	SV_BP_CCF_v3_nohot_tsv\n')

for PAIR_SET=PAIR_SETS
    PAIR_SET
    PS=load_tsv(['FC/' PAIR_SET{1} '.pair_set.tsv']);
    P=load_tsv(['FC/REBC.pairs.' PAIR_FCDATE '.tsv']);
    N=length(PS.pair)
    for i=1:N
        x1=trimStruct(XPASS,ismember(XPASS.individual,PS.pair(i)));
        [PS.pair{i} '   ' num2str(x1.N)]
        file0=['/tmp/SV/' PS.pair{i} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv'];
        file1=['gs://fc-035f5652-acf7-4642-abb7-e8c10848c8ed/SV/' PS.pair{i} '.SV_BP_CCF_v3.noHot.' TODAY '.tsv'];
        printStruct(x1,-1,file0);
        cmd=['gsutil -m cp ' file0 '  ' file1];
        fprintf(fid,'%s %s\n',PS.pair{i},file1);        
    end
end
fclose(fid)

%%
if(0)
    %fragile sites
    %https://webs.iiitd.edu.in/raghava/humcfs/download.html
    
    cmd='grep -E "^chr" SV/FragileSites/fragile_gene.txt > /tmp/fragile_genea.txt'
    unix(cmd)
    cmd='sed ''s/		/	/g''  /tmp/fragile_genea.txt > /tmp/fragile_gene.txt'
    unix(cmd)
    FB=load_bed('/tmp/fragile_gene.txt')
    FB.p1=str2double(FB.p1);
    FB.p2=str2double(FB.p2);
    FB.x1=xhg19(FB.chr,FB.p1);
    FB.x2=xhg19(FB.chr,FB.p2);
    k=find((~isnan(FB.p1))&( ~isnan(FB.p2) ) )
    FB=trimStruct(FB,k);
    [x k]=sort(FB.x1);
    FB=trimStruct(FB,k);
    
    plot(FB.x1,diff([0; FB.x1]))
    
    plot(FB.x2,diff([0; FB.x2]))
    plot(diff([0; FB.x2]))
    plot(FB.x2-FB.x1)
    k=find((FB.x2-FB.x1)>0);
    FB=trimStruct(FB,k);
    k=find((FB.x2-FB.x1)<1e6);
    FB=trimStruct(FB,k);
    semilogy(FB.x2-FB.x1)
    plot(diff([0; FB.x2]))
    
    kx=find(diff([0; FB.x2])<=0);
    k=find(diff([0; FB.x2])>0);
    while length(kx)>0
        length(kx)
        FB=trimStruct(FB,k);
        k=find(diff([0; FB.x2])>0);
        kx=find(diff([0; FB.x2])<=0);
    end
    
    
    XH.x1=xhg19(XH.chr1,XH.pos1);
    XH.x2=xhg19(XH.chr2,XH.pos2);
    
    % good het in TWIST target
    [i,m] = point_overlap_intervals(XH.x1,FB);
    tab(i)
    [i,m] = point_overlap_intervals(XH.x2,FB);
    tab(i)
    
    [i1,m1] = point_overlap_intervals(X.x1,FB);
    tab(i1)
end



