function [G]=ClusterNearestNeighbors1D(x,xc)
% x=rand(100,1); xc=0.01
G=1:length(x);
n=0;
N=length(x);
[x1,m]=sort(x(:));
dx=diff([x1]);
k=(dx>xc);
c=cumsum([1; k]);
G(m)=c;
