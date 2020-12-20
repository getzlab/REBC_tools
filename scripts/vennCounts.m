function [sv] =  vennCounts(X)
[nx ny]=size(X)
N=nx;
b=zeros(N,1);
for j=1:ny
    k=find(X(:,j)>0);
    b(k)=b(k)+2^(j-1);
end
n=hist(b,(1:(2^ny))-1);
sv.data=X;
sv.count=n;
sv.bin=b;
sv.Ncat=ny;
sv.x1=zeros(1,ny);
for d=1:ny
    k=find(sv.data(:,d)==1);sv.x1(d)=sum(sum(sv.data(k,:),2)==1)/length(k);
end
end
function test
f='/Users/stewart/Projects/Cancer/UCS/Fusions/Algs.30Nov2014.txt'
F=load_tsv(f)
sv = vennCounts([F.aMDA F.aUNC F.aBI])
vennDiagram(sv,{'MDA' 'UNC' 'BI'}) 
end
