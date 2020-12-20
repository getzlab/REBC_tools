function tab=load_table(fname,dlm,headlines, getcols, maxlines,varargin)
%function tab=load_table(fname,dlm,headlines, getcols, maxlines, varargin)
%
% returns matlab structure based on delimited tex file
% input:    fname:  file name
%           dlm: optional delimiter (tab = char(9))
%           headlines: default 1 header line with field names
%           getcols: specify indices of fields to return,
%           maxlines: maximum number of lines to load
%
% output:  structure with fields
%

[pathstr, name, ext] = fileparts(fname) ;
matfile= ismember(cellstr(ext),'.mat');
if (matfile)
    if (exist(fname,'file'))
        dm=dir(fname);
        f0=[pathstr '/' name];
        df=dir(f0);
        if isempty(df)
            fprintf('no matching txt file:\t%s\n',f0)
        end
        if (df.datenum<dm.datenum)
            tab=load(fname);
            return;
        else
            fprintf('newer txt file:\t%s\n',f0)
            fprintf('remake mat file:\t%s\n',fname)
            fname = f0;
        end
    else
        fname=fullfile(pathstr,name);
        %fname=[pathstr '/' name];
    end
end


set=exist('dlm','var');
if (set), set=~isempty(dlm); end
if (~set)
    dlm=char(9);
end
set=exist('headlines','var');
if (set), set=~isempty(headlines); end
if (~set)
    headlines=1;
    headlines=autodetectHeader(fname);
    if (headlines>1)
        fprintf('autodetect %d header lines\n',headlines);
    end
end
set=exist('getcols','var');
if (set), set=~isempty(getcols); end
if (~set)
    getcols=[];
end
set_maxlines=exist('maxlines','var');
if (set_maxlines), set_maxlines=~isempty(maxlines); end

q=dir(fname);
if (numel(q)<1)
    error([' no file ' fname]);
end

tab={};
headline={'#'};
fid=fopen(fname,'r');
if (fid<1)
    error([' failed to open ' fname]);
end
% first line
nline=0;
if (headlines<1)
    if (iscell(getcols) )
        % nop header use column names from getcols - must be right
        tab.header=getcols;
        %
        tline = fgetl(fid);
        frewind(fid);
        if ~ischar(tline), return; end
        s=regexp(tline,dlm,'split');
        if (length(s)~=length(getcols))
            fprintf('%s\n',join('\t',getcols));
            fprintf('%s\n',tline);
            error(['getcols do not match columns in file ' fname '\n']);
            return;
        end
        
    else ( (headlines<1)&(~iscell(getcols)) )
        % no header, use v# column names
        tline = fgetl(fid);
        frewind(fid);
        if ~ischar(tline), return; end
        v=regexp(tline,dlm,'split');
        NV=length(v);
        v=cellstr(sprintf('\tv%d',(1:NV))); v=regexp(v,'\t','split'); v=v{1}; v(1)=[];
        tab.header=v;
    end
    headlines=0;
else
    tline=[];
    while (nline<headlines)
        tline = fgetl(fid);
        nline=1+nline;
        headline{nline,1}=tline;
        tab.headline{nline,1}=tline;
    end
    if ~ischar(tline), return; end
    tab.header=regexp(tline,dlm,'split');
    if (~isempty(getcols)&(iscell(getcols)) )
        %tline = fgetl(fid);
        
        [keep,km]=ismember(upper(tab.header),upper(getcols));
        %tab.header=regexp(tline,dlm,'split');
        if (sum(keep)<length(getcols))
            fprintf('%s\n',join('\t',getcols));
            fprintf('%s\n',char(strcat(tab.header))');
            error(['getcols mismatch file: ' fname '\n']);
            return;
        end
        %tab.header=getcols;
        
    end
end
% don't allow blank in last header column
sl=cellfun(@length,tab.header);

if (sl(end)<1)
    tab.header=tab.header(1:(end-1));
    fprintf('%s\n', ['remove tab from end of header line: ' fname ]);
end

NV=length(tab.header);
format=[repmat('%s',1,NV)] ; % \n

fseek(fid,0,'bof');

bufsize=1024*64-1;
cellstrings=false;
if (nargin>5)
    kb=find(ismember(varargin,{'bufsize'}));
    if ~isempty(kb)
        bufsiz=str2int(vargargin{kb+1});
    end
    kb=find(ismember(varargin,{'cellstrings'}));
    cellstrings=~isempty(kb);
end



OK = false;
while ~OK
    try
        if strncmp(char(version),'9.',2)  % matlab 2016a doesn't have bufsize argument...
            tab.dat=textscan(fid,format,'HeaderLines',headlines,'Delimiter',dlm);
        else
            if set_maxlines
                tab.dat=textscan(fid,format,maxlines,'HeaderLines',headlines,'Delimiter',dlm,'BufSize',bufsize); %,varargin{:});
            else
                tab.dat=textscan(fid,format,'HeaderLines',headlines,'Delimiter',dlm,'BufSize',bufsize); %,varargin{:});
            end
        end
        OK = true;
        
    catch exception
        fseek(fid,0,'bof');
        bufsize=(bufsize+1)*2-1;
        if (bufsize>10e7)
            rethrow(exception)
            error(['bufsize: ' num2str(bufsize) '\n']);
            return
        end
        
    end
end
fclose(fid);

NL=min(cellfun(@length,tab.dat((2:end)-1)));
% problem with last line skipping if no end-of-line
% C = textscan(fid, '%s', 'delimiter', sprintf('\n'));C=C{1};C=regexp(C,dlm,'split');
% fix missing CR on last line...
NLend=min(cellfun(@length,tab.dat(end)));
if (NLend==(NL-1))
    tab.dat{end}=[tab.dat{end}; {''}];
end

if (~isempty(getcols)&(~iscell(getcols)) )
    tab.header=tab.header(getcols);
    tab.dat=tab.dat(:,getcols);
    NV=length(tab.header);
end

if (~isempty(getcols)&(iscell(getcols)) )
    [keep,km]=ismember(upper(tab.header),upper(getcols));
    tab.header=tab.header(keep);
    tab.dat=tab.dat(:,keep);
    NV=length(tab.header);
end


for i=1:NV
    v1=tab.dat{i};
    tab.dat{i}=NaN;
    v1=v1(1:NL);
    s1=tab.header{i};
    s1=strtrim(s1);
    s1=regexprep(s1,' ','_');
    s1=regexprep(s1,'(','_');
    s1=regexprep(s1,')','_');
    s1=regexprep(s1,'[','_');
    s1=regexprep(s1,']','_');
    s1=regexprep(s1,'/','_');
    s1=regexprep(s1,'#','');
    s1=regexprep(s1,'%','');
    s1=regexprep(s1,'&','');
    s1=regexprep(s1,'>','_');
    s1=regexprep(s1,'-','_');
    s1=regexprep(s1,'\|','_');
    s1=regexprep(s1,'\.','_');
    s1=regexprep(s1,'__','_');
    s1=regexprep(s1,'__','_');
    s1=regexprep(s1,':','_');
    s1=regexprep(s1,'?','_');
    s1=regexprep(s1,'+','_');
    s1=regexprep(s1,'~','_');
    s1=regexprep(s1,'@','_');
    s1=regexprep(s1,'_$','');
    if (length(s1)<1), continue; end
    if (~isletter(s1(1)))
        s1=['v_' s1];
    end
    
    if (~cellstrings)
        abc=any(cellfun(@length,regexp(v1,'[b-c|f-m|o-z|B-C|F-M|O-Z|:]','start'))>0);
        if (~abc)
            
            if length(which('str2doubleq'))<1
                %x1=str2double(v1);
                x1 = [];
                for a=1:length(v1)
                    if ~isempty(strmatch('NA',v1(a)))
                        x1(a,1) = NaN;
                    else
                        x1(a,1)=str2double(v1(a));
                    end
                end
            else
                vx=v1;
                k=cellfun(@length,v1)<1;
                vx(k)={'NaN'};
                x1=str2doubleq(vx);
            end
            
            k=isnan(x1);
            %abc=length(char(v1(k)))>1;
            %if (sum(k)>0 & abc)
            if (sum(k)>0 )
                abc=~all(ismember(upper(v1(k)),{'NA','NAN','---','\.'}));
            end
        end
        if (~abc)
            v1=x1;
        end
    end
    %     abc=any(cellfun(@length,regexp(regexprep(v1,'NaN','0'),'[b-c|f-m|o-z|B-C|F-M|O-Z/]','start'))>0);
    tab.(s1)=v1;
end
tab=rmfield(tab,'dat');
if (matfile)
    fname_mat=[fname '.mat']
    if ~exist(fname_mat,'file')
        [fid, message]=fopen(fname_mat,'w');
        if fid<0
            fprintf('%s\n',message);
            return;
        else
            fclose(fid)
            save(fname_mat,'-struct','tab');
        end
    end
end


function headlines=autodetectHeader(fname)
commentchars={'#','%','@'};
fid=fopen(fname,'r');
if (fid<1)
    error([' failed to open ' fname]);
end
headlines=1;
tline = fgetl(fid);
while ischar(tline)
    if ~isempty(tline)
        if ~strcmp(tline(1),commentchars)
            fclose(fid);
            return
        end
    end
    headlines=headlines+1;
    tline = fgetl(fid);
end
fclose(fid);



function test
unix('grep -100 Hugo /local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated > ~/test.maf')
X=load_table('~/test.maf',[],[], [1:17 64:65 81:86])



f='/local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated'
f='/local/cga-fh/cga/An_ESO/Individual_Set/PR_Esophageal_CIP_WGS/jobs/wgs/mut/aggregated/PR_Esophageal_CIP_WGS.maf.annotated'


clear
cd ~/Projects/Cancer/DoubleNormal/FH
f1='~/Projects/Cancer/DoubleNormal/FH/LUSC_Samples.tsv'
f1='~/Projects/Cancer/DoubleNormal/FH/test.tsv'
S1=load_table(f1);

f1='/Users/stewart/Projects/Cancer/ExonCapture/refGene.hg19.20100825.sorted.header.txt'
G=load_table(f1,char(9),1);

f2='/Users/stewart/Projects/Cancer/ExonCapture/whole_exome_agilent_1.1_refseq_plus_3_boosters_plus_10bp_padding_minus_mito.Homo_sapiens_assembly19.targets.interval_list.notabtab.txt'
T=load_table(f2,char(9),0,{'chr','p1','p2','strand','targ'});

f3='/Users/stewart/Projects/Cancer/ExonCapture/UCSC.HG19.wgEncodeCrgMapabilityAlign40mer.bwig'
M40=load_table(f3,char(9),1,{'chr','p1','p2','map'});
M40=rmfield(M40,{'dat'});
M40.a=chrom2num(char(M40.chr));



f='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/covbed/EwngSRC-SJDES001.bed'
t=load_table(f,char(9),-1,{'chr','p1','p2','str','targ','nt','nn'});

f='/Users/stewart/Projects/Cancer/Pediatric/Rhabdoid/PR_SIGMA_Rhabdoid.review.clinicalinfo.absolute.maf.txt'
t=load_table(f,char(9));

ln=1;
fl{ln}=dlmsep(tline,dlm);
if ischar(fname)
    fid=fopen(fname,'r');
else
    fid=fname;
end

ln=1;
if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
else
    do_close=0;
end

had_output=0;
fl={};
while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    fl{ln}=dlmsep(tline,dlm);
    ln=ln+1;
    if mod(ln,1000)==0
        verbose(['...' num2str(ln)],30);
        had_output=1;
    end
end
if do_close
    fclose(fid);
    fid=-1;
end
ln=ln-1;
if had_output
    verbose([newline],30);
end



fpos=0;
if headerlines~=0
    if ischar(fname)
        f=fopen(fname,'r');
    else
        f=fname;
        fpos=ftell(f);
    end
    if length(dlm)~=1
        tab.dlm=find_dlm(fname,dlm);
    else
        tab.dlm=dlm;
    end
    if headerlines>0
        tab.headers=read_dlm_file(f,dlm,headerlines);
    elseif headerlines==-1 % R convention
        headerlines=1;
        tab.headers=read_dlm_file(f,dlm,headerlines);
        tab.headers{1}=['EMPTY' tab.headers{1,:}];
    end
    %   fclose(f);
else
    if ischar(fname)
        f=fopen(fname,'r');
    else
        f=fname;
        fpos=ftell(f);
    end
    tab.headers={};
    tab.dlm=dlm;
end

verbose(['Reading file using format:' format],10);
fseek(f,fpos,'bof');
tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
fclose(f);