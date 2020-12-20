function [A, F, A1, F1,AF,DF,segs]=revise_GISTIC_sample_tables(segfile,amp_focal_file,del_focal_file, broad_sig_file, cnv_blacklist_file, arm_coordinates_file,  PAR)
%REVISE_GISTIC_SAMPLE_TABLES generates variant-sample-tables patterned after GISTIC with bugs fixed 
%   Detailed explanation goes here
if nargin<7
    PAR=[];
    PAR.single_amp_threshold = 0.1;
    PAR.single_del_threshold = -0.1;
    PAR.double_amp_threshold = 0.9;
    PAR.double_del_threshold = -0.9;
    PAR.arm_length_fraction_threshold = 0.5;
    PAR.focal_length_fraction_threshold = 0.15;
    PAR.significance_threshold = 0.1; 
    PAR.renormalize=1
end
if nargin<6
    % assume all arms
    arm_coordinates_file=''    
end

if nargin<5
    cnv_blacklist_file=''
end



segs=load_tsv(segfile);
if ~isfield(segs,'Start')
    segs.Start=segs.Start_bp;
    segs.End=segs.End_bp;
    if isfield(segs,'tau')
        segs.Segment_Mean=log2(segs.tau/2);
    else
        segs.Segment_Mean=log2(segs.total_copy_ratio);
    end
    segs.Sample=segs.sample;
    
end

ff=fieldnames(segs);
fkeep={'Sample','Chromosome','Start','End','n_probes','Segment_Mean'};
ffx=ff(~ismember(ff,fkeep));
segs=rmfield(segs,ffx);
segs=orderfields(segs,fkeep);

if PAR.renormalize
    s=unique(sort(segs.Sample));
    segs.CR=2.^(segs.Segment_Mean);
    for i=1:length(s)
        k=find(ismember(segs.Sample,s(i)));        
        meanCR=sum( (segs.End(k)-segs.Start(k)).*segs.CR(k))/sum(segs.End(k)-segs.Start(k));
        segs.CR(k)=segs.CR(k)/meanCR;
        segs.Segment_Mean(k)=log2(segs.CR(k));
    end
end
AF=load_tsv(amp_focal_file);
AF=trimStruct(AF,1:3); AF=rmfield(AF,'N'); AF=flipStruct(AF); N=length(AF.cytoband); AF=trimStruct(AF,1:(N-1));
AF.chr=str2double(cellfun(@(x) x(1),regexp(regexprep(AF.widePeakBoundaries,'chr',''),':','split')));
x=cellfun(@(x) x(2),regexp(AF.widePeakBoundaries,':','split'));
AF.Start=str2double(cellfun(@(x) x(1),regexp(x,'-','split')));
AF.End=str2double(cellfun(@(x) x(2),regexp(x,'-','split')));
AF.gstart=xhg19(AF.chr,AF.Start);
AF.gend=xhg19(AF.chr,AF.End);
AF.cytoband=erase(AF.cytoband,"x");
AF.cytoband=strcat(AF.cytoband,':AMP');
AF.Descriptor=AF.cytoband;


DF=load_tsv(del_focal_file);
DF=trimStruct(DF,1:3); DF=rmfield(DF,'N'); DF=flipStruct(DF); N=length(DF.cytoband); DF=trimStruct(DF,1:(N-1));
DF.chr=str2double(cellfun(@(x) x(1),regexp(regexprep(DF.widePeakBoundaries,'chr',''),':','split')));
x=cellfun(@(x) x(2),regexp(DF.widePeakBoundaries,':','split'));
DF.Start=str2double(cellfun(@(x) x(1),regexp(x,'-','split')));
DF.End=str2double(cellfun(@(x) x(2),regexp(x,'-','split')));
DF.gstart=xhg19(DF.chr,DF.Start);
DF.gend=xhg19(DF.chr,DF.End);
DF.cytoband=erase(DF.cytoband,"x");
DF.cytoband=strcat(DF.cytoband,':DEL');
DF.Descriptor=DF.cytoband;

B=load_tsv(broad_sig_file);
B.sigAMP=zeros(size(B.Arm));
B.sigAMP(find(B.AmpQ_value< PAR.significance_threshold))=1;
B.sigDEL=zeros(size(B.Arm));
B.sigDEL(find(B.DelQ_value< PAR.significance_threshold))=1;


ARMS=load_tsv(arm_coordinates_file);
ARMS.gstart=xhg19(ARMS.chrn,ARMS.start);
ARMS.gend=xhg19(ARMS.chrn,ARMS.xEnd);

%load seg file, trim unnecessarily long names etc.
segs.gstart = xhg19(segs.Chromosome,segs.Start);
segs.gend = xhg19(segs.Chromosome,segs.End);
segs.length = segs.End - segs.Start;

if ~isempty(cnv_blacklist_file)
    %remove segments overlapping with blacklisted regions
    cnv_blacklist = load_tsv(cnv_blacklist_file);
    segs = apply_cnv_blacklist(segs,cnv_blacklist,AF,DF,ARMS);
    segs.length = segs.End - segs.Start;
end

% already renormalized 
%correct segment_means for each sample by sample median
% for samp = unique(segs.Sample)'
%     six = find(ismember(segs.Sample,samp));
%     sampseg = trimStruct(segs,six);
%     sample_median = calc_region_median(sampseg,min(ARMS.gstart),max(ARMS.gend),PAR.arm_length_fraction_threshold);
%     segs.Segment_Mean(six) = segs.Segment_Mean(six) - sample_median;
% end


segs.log_segment_mean = segs.Segment_Mean;
%segs.Segment_Mean = 2.^(segs.Segment_Mean+1)-2;


%generate arm level calls first so we can check focals against them
ARMS=trimStruct(ARMS,ismember(ARMS.GISTICARM,'Y'));

A = [];
A.pair_id=unique(segs.Sample);
A.arm=ARMS.arm';
A.chrn=ARMS.chrn';
A.start=ARMS.start';
A.stop=ARMS.xEnd';
[i m]=ismember(A.arm,B.Arm');
A.sigGAIN(i)=B.sigAMP(m(i));
A.sigDEL(i)=B.sigDEL(m(i));

A1=A;
j=0;
check=0;
for pair=A.pair_id'
    j=j+1;
    segs1=trimStruct(segs,ismember(segs.Sample,pair));
    for i=1:length(ARMS.arm)        
        arm = ARMS.arm{i};
        chrarm = ['chr' arm];
        %select segments that overlap arm, don't select by segment mean yet
        segix = (segs1.gstart <= ARMS.x1(i) & segs1.gend >= ARMS.x1(i)) | (segs1.gstart < ARMS.x2(i) & segs1.gend > ARMS.x2(i)) | (segs1.gstart >= ARMS.x1(i)& segs1.gend <= ARMS.x2(i)  );
        seg1 = trimStruct(segs1,segix);
        %skip event if patient has no supporting segs
        if isempty(seg1.Chromosome) | (sum(seg1.length) < ((ARMS.x2(i) - ARMS.x1(i))*PAR.arm_length_fraction_threshold))
            check=check+1;
            continue
        end
        region_median = calc_region_median(seg1,ARMS.x1(i),ARMS.x2(i),PAR.arm_length_fraction_threshold);
        A1.(chrarm)(j,1) = region_median;
        A.([chrarm '_DEL'])(j,1) = 0;
        if region_median < PAR.single_del_threshold & region_median >= PAR.double_del_threshold
            A.([chrarm '_DEL'])(j,1) = 1;
        elseif region_median < PAR.double_del_threshold
            A.([chrarm '_DEL'])(j,1) = 2;
        end
        
        A.([chrarm '_GAIN'])(j,1) = 0;
        if region_median > PAR.single_amp_threshold & region_median <= PAR.double_amp_threshold;
            A.([chrarm '_GAIN'])(j,1) = 1;
        elseif region_median > PAR.double_amp_threshold;
            A.([chrarm '_GAIN'])(j,1) = 2;
        end
        
        
    end
end
% 

% FOCAL AMPS
F=[];
F.pair_id=unique(segs.Sample);
F.peak = AF.cytoband';
F.peak=regexprep(F.peak,'AMP','GAIN');
F.arm = cellfun(@cellstr,regexp(F.peak,'\d+[p|q]','match'));
F.chr=AF.chr'
F.start=AF.Start';
F.stop=AF.End';

F1=F;
j=0;
check=0;
check3=0;
for pair=F.pair_id'
    j=j+1;
    for alix=1:length(AF.Descriptor)
        peak = AF.Descriptor{alix};
        if ismember(peak(2),{'q';'p'})
            arm = peak(1:2);
        else
            arm = peak(1:3);
        end
        chrarm = ['chr' arm];
        
        %Select segments overlapping with peak
        segix = ((segs.gstart < AF.gstart(alix) & segs.gend > AF.gstart(alix)) | (segs.gstart < AF.gend(alix) & segs.gend > AF.gend(alix)) | (segs.gstart > AF.gstart(alix) & segs.gend < AF.gend(alix))) & ismember(segs.Sample,pair);
        seg1 = trimStruct(segs,segix);
        %only select segments fully enclosed by the peak
        segix2 = (segs.gstart > AF.gstart(alix) & segs.gend < AF.gend(alix)) & ismember(segs.Sample,pair);
        seg2 = trimStruct(segs,segix2);
        
        %length(seg1.gstart)
        if length(seg1.gstart) ==1
            check=check+1;
        end
        if ~isempty(seg1.gstart)
            %Subtract out arm level events ONLY IF PEAK OCCURS IN SIGNIFICANTLY COPY-VARIED ARM
            if (ismember(arm,ARMS.arm) & strfindk(F.peak(alix),'GAIN'))
                armix = find(ismember(ARMS.arm,arm));
                arm_median = A1.(chrarm )(j);
                
                % subtract arm baseline only if arm is significant
                iarm=find(ismember(A.arm,arm)) ;
                arm_median = A1.(chrarm)(j)*A1.sigGAIN(iarm);
                
                seg1.Segment_Mean = seg1.Segment_Mean - arm_median;
            end
            if length(seg1.Chromosome) == 0
                continue
            end
            if strfindk(F.peak(alix),'GAIN')
                seg1.Segment_Mean = -seg1.Segment_Mean;
            end
            region_median = calc_region_median(seg1,AF.gstart(alix),AF.gend(alix),PAR.focal_length_fraction_threshold);
            %correct median back to the sign it originally had
            if strfindk(F.peak(alix),'GAIN')
                region_median = -region_median;
            end
        end
        
        if length(seg1.gstart) == 0
            region_median=NaN; % 'no_coverage';
            j
            pair
            peak
        end
        %create a valid field for the event (amp/del).
        key1 = extractBefore(F.peak{alix},":");
        key = ['chr' regexprep(F.peak{alix},':','_')];
        
        if region_median>-1e9; % ~='no_coverage'
            F1.(key)(j,1) = region_median;
            F.(key)(j,1) = 0;
            if region_median > PAR.single_amp_threshold & region_median <= PAR.double_amp_threshold;
                F.(key)(j,1) = 1;
            elseif region_median > PAR.double_amp_threshold;
                F.(key)(j,1) = 2;
            end
            %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
            if length(seg2.Chromosome) > 0
                segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean > PAR.single_amp_threshold);
                if (sum(segix3) > 0) & (F.(key)(j,1) ==0)
                    F.(key)(j,1) = 1;
                end
            end
        end
        
        if isnan(region_median)
            F1.(key)(j,1) = -999999999999999;
            F.(key)(j,1) = -999999999999999;            
        end
    end
end

FA=F;
FA1=F1;

% 
% FOCAL DELS

F=[];
F.pair_id=unique(segs.Sample);
F.peak = DF.cytoband';
F.arm = cellfun(@cellstr,regexp(F.peak,'\d+[p|q]','match'))
F.chr=DF.chr'
F.start=DF.Start';
F.stop=DF.End';

F1=F;
j=0;
check=0;
check2=0;
check3=0;
for pair=F.pair_id'
    j=j+1;
    for alix=1:length(DF.Descriptor)
        peak = DF.Descriptor{alix};
        if ismember(peak(2),{'q';'p'})
            arm = peak(1:2);
        else
            arm = peak(1:3);
        end
        chrarm = ['chr' arm];
          
        %Select segments overlapping with peak
        segix = ((segs.gstart < DF.gstart(alix) & segs.gend > DF.gstart(alix)) | (segs.gstart < DF.gend(alix) & segs.gend > DF.gend(alix)) | (segs.gstart > DF.gstart(alix) & segs.gend < DF.gend(alix))) & ismember(segs.Sample,pair);
        seg1 = trimStruct(segs,segix);
        %only select segments fully enclosed by the peak
        segix2 = (segs.gstart > DF.gstart(alix) & segs.gend < DF.gend(alix)) & ismember(segs.Sample,pair);
        seg2 = trimStruct(segs,segix2);
        
        if length(seg1.gstart) ==1
            check=check+1;
        end
        if length(seg1.gstart) > 0
            %Subtract out arm level events ONLY IF PEAK OCCURS IN SIGNIFICANTLY COPY-VARIED ARM
            if (ismember(arm,ARMS.arm) & strfindk(F.peak(alix),'DEL'))
                armix = find(ismember(ARMS.arm,arm));
                % subtract arm baseline only if arm is significant
                iarm=find(ismember(A.arm,arm)) ;
                arm_median = A1.(chrarm)(j)*A1.sigDEL(iarm);
                seg1.Segment_Mean = seg1.Segment_Mean - arm_median;
                if seg1.Segment_Mean ==0
                    check2=check2+1;
                end
            end
            if length(seg1.Chromosome) == 0
                continue
            end
            %restrict length only within the peak.
            %subtract those that lay outside of the current peak.
            region_median = calc_region_median(seg1,DF.gstart(alix),DF.gend(alix),PAR.focal_length_fraction_threshold);
            %correct median back to the sign it originally had
            if strfindk(F.peak(alix),'GAIN')
                region_median = -region_median;
            end
        end
        
        if length(seg1.gstart) == 0
            region_median=NaN; %'no_coverage';
            j
            pair
            peak
        end
        %create a valid field for the event (amp/del).
        key1 = extractBefore(F.peak{alix},":");
        key = ['chr' regexprep(F.peak{alix},':','_')];    
        
        if region_median >-1e9
            F1.(key)(j,1) = region_median;
            F.(key)(j,1) = 0;
            if region_median < PAR.single_del_threshold & region_median >= PAR.double_del_threshold
                F.(key)(j,1) = 1;
            elseif region_median < PAR.double_del_threshold
                F.(key)(j,1) = 2;
            end
            %add another or condition (if a segment with (#probe>10, segment_mean>0.1) exists, also set it to be 1)
            if length(seg2.Chromosome) > 0
                segix3 = (seg2.length >= 10000) & (seg2.Segment_Mean < PAR.single_del_threshold);
                if (sum(segix3) > 0) & (F.(key)(j,1) ==0)
                    F.(key)(j,1) = 1;
                end
            end
        end
        
        if isnan(region_median) % =='no_coverage'
            F1.(key)(j,1) = -999999999999999;
            F.(key)(j,1) = -999999999999999;
        end
        
    end
end

fa=fieldnames(FA);
fa=fa(strfindk(fa,'_GAIN'));
for i=1:length(fa)
    F.(fa{i})=FA.(fa{i});
end

F.peak=[F.peak FA.peak];
F.arm=[F.arm FA.arm];
F.chr=[F.chr FA.chr];
F.start=[F.start FA.start];
F.stop=[F.stop FA.stop];
               
f1=fieldnames(F1);
fa1=fieldnames(FA1);
for i=1:length(fa1)
    F1.(fa1{i})=FA1.(fa1{i});
end
F1.peak=[F1.peak FA1.peak];
F1.arm=[F1.arm FA1.arm];
F1.chr=[F1.chr FA1.chr];
F1.start=[F1.start FA1.start];
F1.stop=[F1.stop FA1.stop];


end


%Funciton for calculating segment medians
function region_median = calc_region_median(segs,bound1,bound2,frac)
    [~, ix] = sort(segs.Segment_Mean);
    %sort from lowest to highest on segment mean
    segs = trimStruct(segs,ix);
    %correct length for parts of segment that do not overlap the region of
    %interest
    segs.length = max((min(segs.gend,bound2) - max(segs.gstart,bound1)),0); %+ sum(counts.total_length);
    %find midpoint of all parts of segment that overlap the region, assuming start of region is 0
    segs=trimStruct(segs,segs.length ~=0);
    medianix = (sum(segs.length))*frac;
    pos = 0;
    for i = 1:length(segs.length);
        %must update position before checking if we've passed median
        %also must be outside if statement
        pos = pos+segs.length(i);
        if pos >= medianix;
            region_median = segs.Segment_Mean(i);
            return;
        end
    end
end


function  segs = apply_cnv_blacklist(segs,cnv_blacklist,AF,DF,ARMS)
    cnv_blacklist.gstart = xhg19(cnv_blacklist.Chromosome,cnv_blacklist.Start);
    cnv_blacklist.gend   = xhg19(cnv_blacklist.Chromosome,cnv_blacklist.End  );
    cnv_blacklist.keep_event = zeros(size(cnv_blacklist.Chromosome));
    %arm_level_significance = trimStruct(arm_level_significance,arm_level_significance.significant_amplification | arm_level_significance.significant_deletion);
    AF.length = AF.gend - AF.gstart;
    DF.length = DF.gend - DF.gstart;
    cnv_blacklist.length = cnv_blacklist.gend - cnv_blacklist.gstart;
    ARMS.length = ARMS.gend - ARMS.gstart;
    for i = 1:length(AF.gstart)
        overlapix = find(((cnv_blacklist.gstart <= AF.gstart(i) & cnv_blacklist.gend >= AF.gstart(i)) | (cnv_blacklist.gstart <= AF.gend(i) & cnv_blacklist.gend >= AF.gend(i)) | (cnv_blacklist.gstart >= AF.gstart(i) & cnv_blacklist.gend <= AF.gend(i))) & (cnv_blacklist.length > AF.length(i)*0.01));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
   for i = 1:length(DF.gstart)
        overlapix = find(((cnv_blacklist.gstart <= DF.gstart(i) & cnv_blacklist.gend >= DF.gstart(i)) | (cnv_blacklist.gstart <= DF.gend(i) & cnv_blacklist.gend >= DF.gend(i)) | (cnv_blacklist.gstart >= DF.gstart(i) & cnv_blacklist.gend <= DF.gend(i))) & (cnv_blacklist.length > DF.length(i)*0.01));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
    for i = 1:length(ARMS.gstart)
        overlapix = find(((cnv_blacklist.gstart <= ARMS.gstart(i) & cnv_blacklist.gend >= ARMS.gend(i)) | (cnv_blacklist.gstart <= ARMS.gstart(i) & cnv_blacklist.gend >= ARMS.gend(i)) | (cnv_blacklist.gstart >= ARMS.gstart(i) & cnv_blacklist.gend <= ARMS.gend(i))) & ((cnv_blacklist.length > ARMS.length(i)*0.01)));
        cnv_blacklist.keep_event(overlapix,1) = 1;
    end
    cnv_blacklist = trimStruct(cnv_blacklist,cnv_blacklist.keep_event);
    segs.remove_seg = zeros(size(segs.Chromosome));
    
    for i = 1:length(cnv_blacklist.Chromosome)
        segix = find((segs.gstart <= cnv_blacklist.gstart(i) & segs.gend >= cnv_blacklist.gstart(i)) | (segs.gstart <= cnv_blacklist.gend(i) & segs.gend >= cnv_blacklist.gend(i)) | (segs.gstart >= cnv_blacklist.gstart(i) & segs.gend <= cnv_blacklist.gend(i)) & ~segs.remove_seg) ;
        newsegs = trimStruct(segs,zeros(size(segs.Chromosome)));
        for j = 1:length(segix)
            %if segment spans entire blacklisted region, remove if
            %blacklisted region is >80% of segment length
            %segments
            if (segs.gstart(segix(j)) <= cnv_blacklist.gstart(i) & segs.gend(segix(j)) >= cnv_blacklist.gend(i)) 
                cnv_length = cnv_blacklist.gend(i) - cnv_blacklist.gstart(i);
                seg_length = segs.length(segix(j));
                if (cnv_length/seg_length > 0.8)
                    segs.remove_seg(segix(j)) = 1;
                end
            %if segments spans just start of blacklisted region, set end of
            %segment to start of blacklisted region
            elseif (segs.gstart(segix(j)) <= cnv_blacklist.gstart(i) & segs.gend(segix(j)) >= cnv_blacklist.gstart(i))
                    segs.gend(segix(j))   = cnv_blacklist.gstart(i);
                    segs.End(segix(j)) = cnv_blacklist.Start(i);
            %if segment spans just end of blacklisted region, set start to
            %end of blacklisted region
            elseif (segs.gstart(segix(j))   <= cnv_blacklist.gend(i) & segs.gend(segix(j)) >= cnv_blacklist.gend(i)) 
                    segs.gstart(segix(j))   = cnv_blacklist.gend(i);
                    segs.Start(segix(j)) = cnv_blacklist.End(i);
            %if segment falls entirely within blacklisted region, remove
            %segment entirely
            elseif (segs.gstart(segix(j)) >= cnv_blacklist.gstart(i) & segs.gend(segix(j)) <= cnv_blacklist.gend(i))
                    segs.remove_seg(segix(j)) = 1;
            end
        end
        segs = trimStruct(segs,~segs.remove_seg);
        segs = mergeStruct(newsegs,segs);
    end
end


function dependency_zip
main='revise_GISTIC_sample_tables'
flist = matlab.codetools.requiredFilesAndProducts( main );

mkdir(['/tmp/' main])
cd(['/tmp/' main])
for i=1:length(flist)
    cmd=['scp  "' flist{i} '" .']
    unix(cmd)
end
unix('ls -lath')
end
