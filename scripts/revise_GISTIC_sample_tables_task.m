function  [A, F, A1, F1,AF,DF,segs] = revise_GISTIC_sample_tables_task(segfile,amp_focal_file,del_focal_file, broad_sig_file, cnv_blacklist_file, arm_coordinates_file,ID, ONE_LOG2CR_THRESHOLD,TWO_LOG2CR_THRESHOLD,ARM_LENGTH_FRACTION_THRESHOLD,FOCAL_LENGTH_FRACTION_THRESHOLD, SIGNIFICANCE_THRESHOLD,RENORMALIZE)
%REVISE_GISTIC_SAMPLE_TABLES generates variant-sample-tables patterned after GISTIC with bugs fixed 
%   Detailed explanation goes here
PAR=[]
PAR.single_amp_threshold = ONE_LOG2CR_THRESHOLD;
PAR.single_del_threshold = -1*ONE_LOG2CR_THRESHOLD;
PAR.double_amp_threshold = TWO_LOG2CR_THRESHOLD;
PAR.double_del_threshold = -1*TWO_LOG2CR_THRESHOLD;
PAR.arm_length_fraction_threshold = ARM_LENGTH_FRACTION_THRESHOLD;
PAR.focal_length_fraction_threshold = FOCAL_LENGTH_FRACTION_THRESHOLD;
PAR.significance_threshold = SIGNIFICANCE_THRESHOLD;
PAR.renormalize=RENORMALIZE;
[A, F, A1, F1,AF,DF,segs]=revise_GISTIC_sample_tables(segfile,amp_focal_file,del_focal_file, broad_sig_file, cnv_blacklist_file, arm_coordinates_file,  PAR)
%%
% load('~/Downloads/DLBCL.mat'); cd('/GoogleDrive/Lymphoma/2019/DLBCL/Gene_Sample_Matrix')
% write ARM tables in format broad_values_by_arm.txt and then equivlent
% broad_samples_by_arm with values 0,1,2 for each arm and polarity
B0=load_tsv(broad_values_by_arm_file)
pair_id=A.pair_id
armLAB=fieldnames(A)
B1=flipStruct(rmfield(A1,{'arm','chrn','start','stop','sigGAIN','sigDEL'}))
B1=RenameField(B1,'pair_id','ChromosomeArm')
B1.ChromosomeArm=regexprep(B1.ChromosomeArm,'chr','')
bff=fieldnames(B1); bff=bff(2:end);
for i=1:length(bff)
B1.(bff{i})=2*(2.^B1.(bff{i}))-2;
end
B2=trimStruct(B1,1); B2=rmfield(B2,'N');
B2.ChromosomeArm{1}='ChromosomeArm'
bff=fieldnames(B2)
for i=2:length(bff)
B2.(bff{i})=pair_id(i-1);
end
BV=mergeStruct(B2,B1)
BV
format short
printStruct(BV)

%% Arm 0 1 2's 
B1=flipStruct(rmfield(A,{'arm','chrn','start','stop','sigGAIN','sigDEL'}))
B1=RenameField(B1,'pair_id','ChromosomeArm')
B1.ChromosomeArm=regexprep(B1.ChromosomeArm,'chr','')
bff=fieldnames(B1); bff=bff(2:end);
B2=trimStruct(B1,1); B2=rmfield(B2,'N');
B2.ChromosomeArm{1}='ChromosomeArm'
bff=fieldnames(B2); bff=bff(2:end)
for i=1:length(bff)
  B2.(bff{i})=pair_id(i);
end
% BU header has matlab modified pair_id field names 
% BU first row has original unmodified pair_id names
%    for subsequent unix sed command to delete modified header pair_id names
% BU rows after first row has values
BU=mergeStruct(B2,B1)
printStruct(BU)

%%
% write FOCAL tables in format all_lesions.txt and then equivlent
% focal_samples_by_arm with values 0,1,2 for each arm and polarity
F0=load_tsv(focal_lesions_file)

% original pair_ids
pair_id=F1.pair_id
peak=fieldnames(F1); peak=peak(7:end) 
FT=flipStruct(rmfield(F1,{'peak','chr','start','stop','arm'}))
FT=RenameField(FT,'pair_id','peak')
FT.peak=regexprep(FT.peak,'chr','')
FT.peak=regexprep(FT.peak,'_DEL',':DEL')
FT.peak=regexprep(FT.peak,'_GAIN',':GAIN')
FT.peak=regexprep(FT.peak,'_','.')
bff=fieldnames(FT); bff=bff(2:end);
for i=1:length(bff)
   FT.(bff{i})=2*(2.^FT.(bff{i}))-2;
end
FT1=trimStruct(FT,1); FT1=rmfield(FT1,'N');
FT1.peak{1}='peak'
bff=fieldnames(FT1); bff=bff(2:end)
for i=1:length(bff)
   FT1.(bff{i})=pair_id(i);
end
% FTV header has matlab modified pair_id field names 
% FTV first row has original unmodified pair_id names
%    for subsequent unix sed command to delete modified header pair_id names
% FTV rows after first row has values
FTV=mergeStruct(FT1,FT)
FTV=rmfield(FTV,'N');
FTV
format short
printStruct(FTV)

FT=flipStruct(rmfield(F,{'peak','chr','start','stop','arm'}))
FT=RenameField(FT,'pair_id','peak')
FT.peak=regexprep(FT.peak,'chr','')
FT.peak=regexprep(FT.peak,'_DEL',':DEL')
FT.peak=regexprep(FT.peak,'_GAIN',':GAIN')
FT.peak=regexprep(FT.peak,'_','.')
bff=fieldnames(FT); bff=bff(2:end);
FT1=trimStruct(FT,1); FT1=rmfield(FT1,'N');
FT1.peak{1}='peak'
bff=fieldnames(FT1); bff=bff(2:end)
for i=1:length(bff)
   FT1.(bff{i})=pair_id(i);
end
% FTU header has matlab modified pair_id field names 
% FTU first row afer header has original unmodified pair_id names
%    for subsequent unix sed command to delete modified header pair_id names
% FTU rows after first row has values
FTU=mergeStruct(FT1,FT)
FTU=rmfield(FTU,'N');
FTU
printStruct(FTU)

% original all_lesions.txt file for format of first few rows... not pair_id fields 
FF=F0
xpairs=fieldnames(FF); xpairs=xpairs(10:end)
FF=rmfield(FF,xpairs)
k=strfindK(FF.AmplitudeThreshold,'^0:')
FFU=trimStruct(FF,k)
FFU=rmfield(FFU,'N');
k=strfindK(FF.AmplitudeThreshold,'^Actual')
FFV=trimStruct(FF,k)
FFV=rmfield(FFV,'N');

pairs=fieldnames(FTU); pairs=pairs(2:end)
FFU.peak=regexprep(strcat(FFU.Descriptor,':DEL'),' ','')
k=strfindK(FFU.UniqueName,'^Amplification')
FFU.peak(k)=regexprep(FFU.peak(k),':DEL',':GAIN')
[i1 m1]=ismember(FFU.peak,FTU.peak)
for i=1:length(pairs)
    %FFU.(pairs{i})=zeros(size(FFU.qValues));
    FFU.(pairs{i})(1) =  FTU.(pairs{i})(1);
    FFU.(pairs{i})(i1,1) =  FTU.(pairs{i})(m1(i1));
end
[i1 m1]=ismember(FTU.peak,FFU.peak)
FTU.UniqueName(i1,1)=regexprep(FFU.UniqueName(m1(i1)),' ','');
FTU.Descriptor(i1,1)=regexprep(FFU.Descriptor(m1(i1)),' ','');
FTU.WidePeakLimits(i1,1)=regexprep(FFU.WidePeakLimits(m1(i1)),' ','');
FTU.PeakLimits(i1,1)=regexprep(FFU.PeakLimits(m1(i1)),' ','');
FTU.RegionLimits(i1,1)=regexprep(FFU.RegionLimits(m1(i1)),' ','');
FTU.qValues(i1,1)=regexprep(cellstr(num2str(FFU.qValues(m1(i1)))),' ','');
FTU.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks(i1,1)=regexprep(cellstr(num2str(FFU.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks(m1(i1)))),' ','');
FTU.BroadOrFocal(i1,1)=regexprep(FFU.BroadOrFocal(m1(i1)),' ','');
FTU.AmplitudeThreshold(i1,1)=regexprep(FFU.AmplitudeThreshold(m1(i1)),' ','');

FTU.UniqueName{1}='UniqueName'
FTU.Descriptor{1}='Descriptor'
FTU.WidePeakLimits{1}='WidePeakLimits'
FTU.PeakLimits{1}='PeakLimits'
FTU.RegionLimits{1}='RegionLimits'
FTU.qValues{1}='qValues'
FTU.UniqueName{1}='UniqueName'
FTU.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks{1}='ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks'
FTU.BroadOrFocal{1}='BroadOrFocal'
FTU.AmplitudeThreshold{1}='AmplitudeThreshold'

FUF=orderfields(FTU,FFU)
peak=FUF.peak;
FUF=rmfield(FUF,'peak')
%FUF.peak=peak;
printStruct(FUF)

pairs=fieldnames(FTV); pairs=pairs(2:end)
FFV.peak=regexprep(strcat(FFV.Descriptor,':DEL'),' ','')
k=strfindK(FFV.UniqueName,'^Amplification')
FFV.peak(k)=regexprep(FFV.peak(k),':DEL',':GAIN')
[i1 m1]=ismember(FFV.peak,FTV.peak)
for i=1:length(pairs)
    FFV.(pairs{i})(1) =  FTV.(pairs{i})(1);
    FFV.(pairs{i})(i1,1) =  FTV.(pairs{i})(m1(i1));
end
[i1 m1]=ismember(FTV.peak,FFV.peak)
FTV.UniqueName(i1,1)=regexprep(FFV.UniqueName(m1(i1)),' ','');
FTV.Descriptor(i1,1)=regexprep(FFV.Descriptor(m1(i1)),' ','');
FTV.WidePeakLimits(i1,1)=regexprep(FFV.WidePeakLimits(m1(i1)),' ','');
FTV.PeakLimits(i1,1)=regexprep(FFV.PeakLimits(m1(i1)),' ','');
FTV.RegionLimits(i1,1)=regexprep(FFV.RegionLimits(m1(i1)),' ','');
FTV.qValues(i1,1)=regexprep(cellstr(num2str(FFV.qValues(m1(i1)))),' ','');
FTV.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks(i1,1)=regexprep(cellstr(num2str(FFV.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks(m1(i1)))),' ','');
FTV.BroadOrFocal(i1,1)=regexprep(FFV.BroadOrFocal(m1(i1)),' ','');
FTV.AmplitudeThreshold(i1,1)=regexprep(FFV.AmplitudeThreshold(m1(i1)),' ','');

FTV.UniqueName{1}='UniqueName'
FTV.Descriptor{1}='Descriptor'
FTV.WidePeakLimits{1}='WidePeakLimits'
FTV.PeakLimits{1}='PeakLimits'
FTV.RegionLimits{1}='RegionLimits'
FTV.qValues{1}='qValues'
FTV.UniqueName{1}='UniqueName'
FTV.ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks{1}='ResidualQValuesAfterRemovingSegmentsSharedWithHigherPeaks'
FTV.BroadOrFocal{1}='BroadOrFocal'
FTV.AmplitudeThreshold{1}='AmplitudeThreshold'

FVF=orderfields(FTV,FFV)
peak=FVF.peak;
FVF=rmfield(FVF,'peak')
%FUF.peak=peak;
printStruct(FVF)


Fv=trimStruct(FVF,(2:length(FVF.Descriptor)))
FL=mergeStruct(FUF,Fv)

% write output files 

BVname=[ID '.broad_values_by_arm.V3.tsv']
printStruct(BV,-1,'/tmp/BV.tsv')
cmd=['sed ''1d'' /tmp/BV.tsv  > ' BVname]
unix(cmd)

BUname=[ID '.broad_bin_by_arm.V3.tsv']
printStruct(BU,-1,'/tmp/BU.tsv')
cmd=['sed ''1d'' /tmp/BU.tsv  > ' BUname]
unix(cmd)

FLname=[ID '.focal_lesions.V3.tsv']
printStruct(FL,-1,'/tmp/FL.tsv')
cmd=['sed ''1d'' /tmp/FL.tsv  > ' FLname]
unix(cmd)

end


function test
%%
clear
cd('/GoogleDrive/Lymphoma/2019/DLBCL/Gene_Sample_Matrix')
%cohort level seg file
ID='DLBCL_304'
segfile = 'sampledata/inputs_for_gene_sample_matrix/DLBCL.allelic_capseg.seg.txt';
amp_focal_file = 'sampledata/DLBCL.report_gistic2/amp_genes.txt'
del_focal_file = 'sampledata/DLBCL.report_gistic2/del_genes.txt'
broad_sig_file = 'sampledata/DLBCL.report_gistic2/broad_significance_results.txt'
broad_values_by_arm_file = 'sampledata/DLBCL.report_gistic2/broad_values_by_arm.txt'
focal_lesions_file = 'sampledata/DLBCL.report_gistic2/all_lesions.txt'

cnv_blacklist_file = 'sampledata/inputs_for_gene_sample_matrix/CNV.hg19.bypos.111213.CR1_event_added.bed';

% arms-coordinates file
arm_coordinates_file   = 'hg19.GISTIC.arms.tsv';
arm_coordinates_file   = 'sampledata/hg19.GISTIC.arms.tsv';
addpath('/GoogleDrive/Cancer/tools/matlab')

ONE_LOG2CR_THRESHOLD = 0.1
TWO_LOG2CR_THRESHOLD = 0.9
ARM_LENGTH_FRACTION_THRESHOLD = 0.5
FOCAL_LENGTH_FRACTION_THRESHOLD =0.15
SIGNIFICANCE_THRESHOLD = 0.1
RENORMALIZE = 1

[A, F, A1, F1, AF, DF]=revise_GISTIC_sample_tables_task(segfile,amp_focal_file,del_focal_file,broad_sig_file,cnv_blacklist_file,arm_coordinates_file,broad_values_by_arm_file,focal_lesions_file,ID, ONE_LOG2CR_THRESHOLD,TWO_LOG2CR_THRESHOLD,ARM_LENGTH_FRACTION_THRESHOLD,FOCAL_LENGTH_FRACTION_THRESHOLD, SIGNIFICANCE_THRESHOLD,RENORMALIZE)

% compare ARMs 

B0=load_tsv(broad_values_by_arm_file)
B1=load_tsv([ID '.broad_values_by_arm.V3.tsv'])

[B0.ChromosomeArm B1.ChromosomeArm]
B1=orderfields(B1,B0)
pairs=fieldnames(B1); pairs=pairs(2:end)
x0=[]
x1=[];
p=cellfun(@(x) x(1),regexp(pairs,'_TP_','split'));
p=regexprep(p,'_nullpair','');
p=regexprep(p,'^DLBCL_','');
p=regexprep(p,'^DFCI_','');
p=regexprep(p,'^DLBCL','')    

figure(1)

for i=1:length(pairs)
    
    X0=B0.(pairs{i});
    X1=B1.(pairs{i})
    DX=X1-X0
    plot(X0,X1,'o')
    k=find(abs(DX)>0.5)
    text(X0(k)+0.05,X1(k)+0.05,strcat(B0.ChromosomeArm(k),':',p(i)))
    x0=[x0;X0] ;
    x1=[x1;X1] ;
    hold on
end
hold off
grid on
xlabel('GISTIC CR*2-2')
ylabel('FIX CR*2-2')


figure(2)
myhist(x0-x1,100,'log')

%% focals


FA0=load_tsv(amp_focal_file)
B1=load_tsv([ID '.broad_values_by_arm.V3.tsv'])

[B0.ChromosomeArm B1.ChromosomeArm]
B1=orderfields(B1,B0)
pairs=fieldnames(B1); pairs=pairs(2:end)
x0=[]
x1=[];
p=cellfun(@(x) x(1),regexp(pairs,'_TP_','split'));
p=regexprep(p,'_nullpair','');
p=regexprep(p,'^DLBCL_','');
p=regexprep(p,'^DFCI_','');
p=regexprep(p,'^DLBCL','')    

figure(1)

for i=1:length(pairs)
    
    X0=B0.(pairs{i});
    X1=B1.(pairs{i})
    DX=X1-X0
    plot(X0,X1,'o')
    k=find(abs(DX)>0.5)
    text(X0(k)+0.05,X1(k)+0.05,strcat(B0.ChromosomeArm(k),':',p(i)))
    x0=[x0;X0] ;
    x1=[x1;X1] ;
    hold on
end
hold off
grid on
xlabel('GISTIC CR*2-2')
ylabel('FIX CR*2-2')



end

%%
function dependency_zip
main='revise_GISTIC_sample_tables_task'
flist = matlab.codetools.requiredFilesAndProducts( main );

mkdir(['/tmp/' main])
cd(['/tmp/' main])
for i=1:length(flist)
    cmd=['scp  "' flist{i} '" .']
    unix(cmd)
end
unix('ls -lath')
end


