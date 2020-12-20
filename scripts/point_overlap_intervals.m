function [i,m] = point_overlap_intervals(a,B)
% point_overlap_intervals returns subset of A that intersects B
% based on a and B.x1, B.x2 values
% x1 x2 must be sorted (no overlaps in B intervals)
% i is 0,1 map of A that overlap with B
% m is index of B that overlaps A in the convention of ismember
i=zeros(size(a));
m=i;

% a after B x1
[n1 b1] = histc(a,[B.x1; B.x1(end)+1]);
% missing B's
% B=trimStruct(B,n11>0)
%[n11 b11] = histc(A.x1,B.x1);
% a after B x2
[n2 b2] = histc(a,[B.x2; B.x2(end)+1]);
k=find(b1==b2+1);
i(k)=1;
i=i>0;
m=0*i;
m(k)=b1(k);
end

