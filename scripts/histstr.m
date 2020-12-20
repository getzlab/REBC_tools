function [ n,b ] = histstr(A,b0,w)
% function [ n,b ] = histstr(A,b0,w)
% histogram occurances of cell strings 
% inputs: 
%   A is list of cell strings
%   b0 list of unique cell string categories (optional) 
%   w vector of weights of same dimension as A (optional) 
% output: 
%   n counts of occurances of b0 in A (weighted counts if w input)
% no output:
%   frequency bar plot
% 
% CS 2016
%
if isnumeric(A)
    A=num2str(A);
end
if (nargin<2)   
    b=sort(unique(A));
else
    b=b0;
end
b=b(:);
n=countmember(b,A);
if exist('w','var')
    [i m]=ismember(A,b);
    um=unique(m);
    nw=zeros(size(b));
    %ii=1:length(A);
    
    for m1=um(:)'
        nw(m1)=nansum(w(m==m1));
    end
    n=nw;
end
if nargout<1
    if length(n)<10
       bar(n);
       set(gca,'xtick',1:length(n),'XTickLabel',b)
     else
       barh(n)
       set(gca,'ytick',1:length(n),'YTickLabel',b)
    end
    set(gca,'TickLabelInterpreter','none')
     
end

end

function test
    a={'Missense_Mutation'
    'Missense_Mutation'
    'Non-coding_Transcript'
    'Silent'
    'Missense_Mutation'
    'IGR'
    'Silent'
    'Missense_Mutation'
    'Intron'
    'Silent'
    'Missense_Mutation'
    'Missense_Mutation'
    'Missense_Mutation'
    'Non-coding_Transcript'
    'Missense_Mutation'
    'Missense_Mutation'
    }
    w=1:16;
    
histstr(a)
b=unique(a)
histstr(a,b,w)
[n b]=histstr(a)
[n b]=histstr(a)
b1=[b;{'fdafsd'}]
[n]=histstr(a,b1)
histstr(a,b1)
b1=[b;{'Aafsd'}]
histstr(a,sort(b1))

end