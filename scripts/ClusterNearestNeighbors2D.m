function [G]=ClusterNearestNeighbors2D(xy,xyc)
%   ClusterNearestNeighbors2D(xy,xyc)  
%   cluster in two dimensions of xy: x(:,1) 
%   with neighbor threshold xyc(1), and  xy(:,2) with threshold xyc(2)
%   inputs:
%       xy = xy(:,1) and xy(:,2) variables to be clustered
%       xyc = xc(1) and xc(2) cluster scale neighbor threshold 
%   outputs: 
%       G = group id for each cluster in order of x
%
% Author: Chip Stewart
% Created: 2007-12-04
% Copyright 2007 Boston College.
%
% x=rand(100,2); x1c=0.01;  x2c=0.02

G=[];
if isempty(xy), return; end
[N,dim]=size(xy);
if (dim~=2), return; end
G=1;
if (N<2), return; end
[Nc]=prod(size(xyc));
if (Nc>2), return; end
if (Nc==2), xs=xyc(1); ys=xyc(2); end 
if (Nc==1), xs=xyc(1); ys=xs;  end
if (Nc<1),  return; end

%cluster first dimension
%x1c=xc(1)
[x1,m]=sort(xy(:,1));
dx=diff([x1]);
k=(dx>=xs);
c=cumsum([1; k]);
G1(m)=c;

% cluster second dimension
%x2c=xc(2)
[y1,m]=sort(xy(:,2));
dy=diff([y1]);
k=(dy>=ys);
c=cumsum([1; k]);
G2(m)=c;

% combine groups
Gx = 10^(1+ceil(log10(max(G1))));
G0 = G2*Gx+G1;

% remap to consecutive group id's 
UG=unique(G0);
G=0*x1;
for n=1:length(UG)
    k=find(G0==UG(n));
    G(k)=n+0*k;
end

if (0)
    UG=unique(G);
    for i=1:length(UG)
        fprintf(1,'%d ',i)
        k=find(G==UG(i));
        fprintf(1,'%.3f ',G(k))
        fprintf(1,'\n')
    end
end
