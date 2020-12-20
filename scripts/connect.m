function [class,nxt]=connect(i,j,class,nxt)
% Saul & Knuth Equivalence Class algorithm
if (class(i)==class(j)), return;end
j1=j;
class(j1)=class(i);
while (nxt(j1)~=j)
    j1=nxt(j1);
    class(j1)=class(i);
end
i1=nxt(i);
nxt(i)=j;
nxt(j1)=i1;
