function [T]=tab(x,opt)
    Q=tabulate(x);
    T.x=Q(:,1);
    T.n=Q(:,2);
    if ~isnumeric(T.n)
        T.n=cell2mat(T.n);
    end
    if nargin<2
        opt='';
    end
    if (~isempty(strfind(upper(opt),'S'))|(nargout<1)|isnumeric(opt))
            [n k]=sort(T.n,1,'ascend');        
            T=trimStruct(T,k);
    end
    if ~isempty(strfind(upper(opt),'R'))
            [n k]=sort(T.n,1,'descend');        
            T=trimStruct(T,k);
    end
    if nargout<1
        if isnumeric(opt)
            n=opt(1);
        else
            n=str2double(regexp(opt,'\d+','match'));
        end
        if isempty(n)
            n=length(T.x);
        end
        printStruct(T,1:n)
        fprintf('---\ntot\t%d\n',sum(T.n))
        fprintf('case\t%d\n',length(T.n))
        clear T
    end
end
