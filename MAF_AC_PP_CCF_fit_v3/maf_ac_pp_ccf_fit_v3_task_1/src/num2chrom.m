function [chr]=num2chrom(a)
%  [chr]=num2chrom(a)
% convert numbered chromosomes to equivalent character names 
% for human assemblies up to 25=MT. Other than X,Y,MT, the number is simply
% converted to ascii. 
% Chip Stewart, 27 May 2011 
%
chr=cell(size(a));
for i=1:numel(a)
    chr{i}=num2str(a(i));
    if (a(i)==23), chr{i}='X'; end;
    if (a(i)==24), chr{i}='Y'; end;
    if (a(i)==25), chr{i}='MT'; end;
end
