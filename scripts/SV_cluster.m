function [X2,Z] = SV_cluster(X,WB,WC)
%[X] = SV_cluster(X,W)
%   cluster SV events in BP format within window W
%   added fields:
%   balanced>0,   SV is part of balanced event pair (two alleles), bal index of pairs
%   inversion2>0,  SV is part of complete inversion break (FF,RR) at boths end of inversion
%   complex>0,  SV is part of complex set of SVs, complex is index of sets
%

if nargin<2
    % window bp for clustered breaks
    WB=1000
    WC=25000
end

X.ID=X.individual;
% if isfield(X,'LABELED_REBC_ID')
%     X.ID=X.LABELED_REBC_ID;
% end
tab(X.ID)
ID=sort(unique(X.ID))
n=histstr((X.ID),ID)
IDX=ID(n>1) % 183 

Z=[]
Z.ID=ID;
X1=[]
X.balanced=0*X.pos1;
X.inversion2=0*X.pos1;
X.complex=0*X.pos1;

if ~isfield(X,'x1')
    X.x1=xhg19(X.chr1,X.pos1);
    X.x2=xhg19(X.chr2,X.pos2);
end

W=WB;
% balanced SVs and inversions with two sets of breakpoints
for i=1:length(ID)
    x1=trimStruct(X,ismember(X.ID,Z.ID(i)));
    x1.clusterDeltaPos1=NaN*x1.pos1;
    x1.clusterDeltaPos2=NaN*x1.pos1;
    Z.nSV(i,1)=x1.N;
    Z.nSVdeletion(i,1)=sum(ismember(x1.class,'deletion'));
    Z.nSVtandem_dup(i,1)=sum(ismember(x1.class,'tandem_dup'));
    Z.nSVinversion(i,1)=sum(ismember(x1.class,'inversion'));
    Z.nSVlong_range(i,1)=sum(ismember(x1.class,'long_range'));
    Z.nSVinter_chr(i,1)=sum(ismember(x1.class,'inter_chr'));
    %if x1.N<2, continue; end
    % look for balanced events
    %bal=ismember(x1.balanced,'True');
    xy=[x1.x1 x1.x2];
    [G]=ClusterNearestNeighbors2D(xy,W*[1 1]);
    t=tab(G);
    k=find(ismember(G,t.x(t.n>1)));
    %[x1.chr1(k) x1.pos1(k) x1.str1(k) x1.chr2(k) x1.pos2(k) x1.str2(k) G(k) ]
    G1=unique(G(k));
    nb=0;
    ni=0;
    for g=1:length(G1)
        k=find(ismember(G,G1(g)));
        kbal=k(find((x1.str1(k)~=x1.str2(k))|(x1.chr1(k)~=x1.chr2(k))));
        if length(kbal)>1
            k1=kbal;            
            bal=find(x1.str1(k1(1))==(1-x1.str1(k1(2))))&(x1.str2(k1(1))==(1-x1.str2(k1(2))));
            if bal
                if length(k1)>2
                    % extra balanced SV ... take out redundant SV with lowest (NALG*100+TALT)
                    tx=tab(x1.str1(k1))
                    sx=tx.x(tx.n>1)
                    kx=k1(find(x1.str1(k1)==sx))
                    [~,kx1]=min(x1.NALG(kx)*100+x1.VCF_TALT(kx))
                    k1(kx(kx1)==k1)=[];
                end
                nb=nb+1;
                x1.balanced(k1)=nb;
                % convention pos1 pos2 check
                if length(unique(x1.pos1(k1)<x1.pos2(k1)))>1
                    print('oops')
                end
                % scar length rev-fwD
                kF=k1(find(x1.str1(k1)==0));  kR=k1(find(x1.str1(k1)==1));

                x1.clusterDeltaPos1(k1)=x1.pos1(kF)-x1.pos1(kR);
                kF=k1(find(x1.str2(k1)==0));  kR=k1(find(x1.str2(k1)==1));
                x1.clusterDeltaPos2(k1)=x1.pos2(kF)-x1.pos2(kR);
            else
                ID{i}, [x1.chr1(k1) x1.pos1(k1) x1.str1(k1) x1.chr2(k1) x1.pos2(k1) x1.str2(k1)]
                %keyboard
            end
        end
        kinv=k(find((x1.str1(k)==x1.str2(k))&(x1.chr1(k)==x1.chr2(k))));
        if length(kinv)>1
            k1=kinv;
            inv=find(x1.str1(k1(1))==(1-x1.str1(k1(2))))&(x1.str2(k1(1))==(1-x1.str2(k1(2))));
            if inv
               if length(k1)>2
                    % extra inversion2 SV ... take out redundant SV with lowest (NALG*100+TALT)
                    tx=tab(x1.str1(k1))
                    sx=tx.x(tx.n>1)
                    kx=k1(find(x1.str1(k1)==sx))
                    [~,kx1]=min(x1.NALG(kx)*100+x1.VCF_TALT(kx))
                    k1(kx(kx1)==k1)=[];
                end                
                ni=ni+1;
                x1.inversion2(k1)=ni;
                kF=k1(find(x1.str1(k1)==0));  kR=k1(find(x1.str1(k1)==1));
                x1.clusterDeltaPos1(k1)=x1.pos1(kF)-x1.pos1(kR);
                kF=k1(find(x1.str2(k1)==0));  kR=k1(find(x1.str2(k1)==1));
                x1.clusterDeltaPos2(k1)=x1.pos2(kF)-x1.pos2(kR);
                % convention pos1 pos2 check
                if length(unique(x1.pos1(k1)<x1.pos2(k1)))>1
                    print('oops')
                end                
                
            else
                ID{i}, [x1.chr1(k1) x1.pos1(k1) x1.str1(k) x1.chr2(k1) x1.pos2(k1) x1.str2(k1)]
                %keyboard
            end
        end
    end
    Z.nSVbal(i,1)=nb;
    Z.nSVinv2(i,1)=ni;
    if isempty(X1)
        X1=x1;
    else
        X1=mergeStruct(X1,x1);
    end
        
end

X2=[];
W=WC;
% chains of SVs 
for i=1:length(ID)
    ID{i}
    x1=trimStruct(X1,ismember(X1.ID,Z.ID(i)));
    class=1:x1.N;
    nxt=class;  % self-linked list
    % collapse balanced and inversion2's
    %G1=ClusterNearestNeighbors1D(x1.x1,W);
    %G2=ClusterNearestNeighbors1D(x1.x2,W);
    
    x=[x1.x1; x1.x2];
    j=[1:x1.N 1:x1.N]';
    G=ClusterNearestNeighbors1D(x,W);

    for i0=1:length(x)
        i1=j(i0);
        k0=find(G==G(i0));
        j1=j(k0)';
        j1(j1==i1)=[];
        for i2=j1   
            [i1 i2];
            [class,nxt]=connect(i1,i2,class,nxt);
        end
    end
    
    GU=unique(class);
    NG=length(GU);
    G=0*class;
    nc=0;
    nSVc=0;
    nSVbal1=sum(x1.balanced>0);
    nSVinv1=sum(x1.inversion2>0);
    for i1=1:NG
        G(class==GU(i1))=i1;
        k1=find(G==i1);
        if length(k1)>1
            if (length(k1)==2)&any(x1.balanced(k1)>0)
                if x1.balanced(k1(1))==x1.balanced(k1(2))
                    continue
                end
            end
            if (length(k1)==2)&any(x1.inversion2(k1)>0)
                if x1.inversion2(k1(1))==x1.inversion2(k1(2))
                    continue
                end
            end
            [x1.chr1(k1) x1.pos1(k1) x1.str1(k1) x1.chr2(k1) x1.pos2(k1) x1.str2(k1) x1.balanced(k1) x1.inversion2(k1)]         
            nc=nc+1;
            nSVc=nSVc+length(k1);
            x1.complex(k1)=nc;
            
            for k=k1(:)'
                if x1.balanced(k), continue; end
                if x1.inversion2(k), continue; end
                
                kx=k1; kx(kx==k)=[];
                [d1x,k1x]=min(abs(x1.x1(k)-x1.x1(kx)));
                [d2x,k2x]=min(abs(x1.x2(k)-x1.x2(kx)));                
                [d12x,k12x]=min(abs(x1.x1(k)-x1.x2(kx)));                
                [d21x,k21x]=min(abs(x1.x2(k)-x1.x1(kx)));    
                d1x=min([d1x d12x])
                d2x=min([d2x d21x])
                if (d1x<W)
                    x1.clusterDeltaPos1(k)=d1x;
                end
                if (d2x<W)
                    x1.clusterDeltaPos2(k)=d2x;
                end
            end
            kb=k1(find(x1.balanced(k1)));
            ki=k1(find(x1.inversion2(k1)));
            nSVbal1=nSVbal1-length(kb);
            nSVinv1=nSVinv1-length(ki);            
        end
    end
    if isempty(X2)
        X2=x1;
    else
        X2=mergeStruct(X2,x1);
    end
    Z.nComplex(i,1)=nc;
    Z.nSVcomplex(i,1)=nSVc;
    Z.nSVbal_notComplex(i,1)=nSVbal1;
    Z.nSVinv2_notComplex(i,1)=nSVinv1;
    
end

printStruct(Z)


end

function test

clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
PSET=PAIR_SETS{1}
unix(['ls -latr SV/' PSET '*.tsv'])

X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filt.24May2020.tsv'])
X0=load_tsv(['SV/' PSET '.aggregated.all.clustered.26May2020.tsv'])

WB=1000
WC=25000
[X1,Z] = SV_cluster(X0,WB,WC)

end


function test1

clear
cd('/GoogleDrive/Cancer/REBC/')
PAIR_SETS={'REBC_primary_pair_393','REBC-NT-NB_240','THCA-PRIMARY-7Dec2017'};
PSET=PAIR_SETS{1}
unix(['ls -latr SV/' PSET '*.tsv'])

X=load_tsv(['SV/' PSET '.SV_BP_CCF_v3.noHot.filt.24May2020.tsv'])
%X0=load_tsv(['SV/' PSET '.aggregated.all.clustered.26May2020.tsv'])

X0=trimStruct(X,ismember(X.LABELED_REBC_ID,'REBC-ADKY'))
X0=printStruct(X0,(X0.chr1==3)|(X0.chr1==10))

WB=1000
WC=25000
[X1,Z] = SV_cluster(X0,WB,WC)

printStruct(X1,~X1.complex)

end


