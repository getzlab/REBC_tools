function X=fc_SV_CCF_v3(X1,ST1,gender,id)

if ischar(X1)
    X1=load_tsv(X1);
end
if ischar(ST1)
    ST1=load_tsv(ST1);
end

fprintf('id: %s',id)
fprintf('purity: %.3f',ST1.purity(1))
fprintf('#SVs: %d',length(X1.pos1))


X=SV_CCF_v3(X1,ST1,gender);


printStruct(X,-1,[id '.SV_CCF_v3.tsv']);

if isdeployed
    exit
end

end

%%
function make1
%%
main='fc_SV_CCF_v3'
flist = matlab.codetools.requiredFilesAndProducts( main );

mkdir(['/tmp/' main])
cd(['/tmp/' main])
for i=1:length(flist)
    cmd=['scp  "' flist{i} '" .']
    unix(cmd)
end
unix('ls -lath')
end