function [B]=load_bed(fname)
%B1=load_table(fname,char(9),0,[],1)


B1 = readtable(fname,'Delimiter',char(9),'FileType','text','HeaderLines',1,'ReadVariableNames',0,'Datetime','text'); %,'Format','auto');
B1 = table2struct(B1,'ToScalar',true);

B.chr=B1.Var1;
B.p1=B1.Var2;
B.p2=B1.Var3;
if isfield(B1,'Var4')
    B.name=B1.Var4;
end
if isfield(B1,'Var5')
    B.score=B1.Var5;
end




function test
Q=load_bed('hg19new.bed');