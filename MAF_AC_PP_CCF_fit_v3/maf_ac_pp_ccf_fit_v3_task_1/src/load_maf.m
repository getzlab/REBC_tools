function X = load_maf(maf_file,OPT)
if nargin>1
    if strfind(OPT,'trim')
        cmd=[' grep Hugo_Symbol ' maf_file];
        [o r]=unix(cmd);
        if length(r)>1
            r=regexp(r,'\t','split');
            ft={
                'Hugo_Symbol'
                'Entrez_Gene_Id'
                'Center'
                'NCBI_Build'
                'Chromosome'
                'Start_position'
                'End_position'
                'Strand'
                'Variant_Classification'
                'Variant_Type'
                'Reference_Allele'
                'Tumor_Seq_Allele1'
                'Tumor_Seq_Allele2'
                'dbSNP_RS' 
                'dbSNP_Val_Status'
                'Tumor_Sample_Barcode'
                'Matched_Norm_Sample_Barcode'
                'Match_Norm_Seq_Allele1'
                'Match_Norm_Seq_Allele2'
                'Tumor_Sample_UUID'
                'Matched_Norm_Sample_UUID'
                'Protein_Change'
                'ref_context'
                'gc_content'
                'i_COSMIC_n_overlapping_mutations'
                'i_dbNSFP_1000Gp1_AF'
                'i_t_lod_fstar'
                'i_tumor_f'
                't_alt_count'
                't_ref_count'
                'i_n_alt_count'
                'i_n_alt_count_realign'
                'i_n_ref_count'
                'i_n_ref_count_realign'
                'i_t_alt_count_realign'
                'i_t_ref_count_realign'
                'i_judgement'
                'alt'
                'ref'
                'q_hat'
                'HS_q_hat_1'
                'HS_q_hat_2'
                'Pr_somatic_clonal'
                'clonal_ix'
                'clonal.ix'
                'homozygous.ix'
                'homozygous_ix'
                'cancer_cell_frac'
                'ccf_CI95_low'
                'ccf_CI95_high'
                'purity'
                'ploidy'
                'i_t_lod_fstar'
                'i_CENTER'
                'i_BCGSC'
                'i_t_alt_count_BCGSC'
                'i_t_ref_count_BCGSC'
                'i_t_depth_BCGSC'
                'i_n_depth_BCGSC'
                'i_n_alt_countBCGSC'
                'i_BCM'
                'i_t_alt_count_BCM'
                'i_t_ref_count_BCM'
                'i_t_depth_BCM'
                'i_n_depth_BCM'
                'i_n_alt_countBCM'
                'i_BI'
                'i_t_alt_count_BI'
                'i_t_ref_count_BI'
                'i_t_depth_BI'
                'i_n_depth_BI'
                'i_n_alt_countBI'
                'i_WU'
                'i_t_alt_count_WU'
                'i_t_ref_count_WU'
                'i_t_depth_WU'
                'i_n_depth_WU'
                'i_n_alt_countWU'
                'i_UCSC'
                'i_t_alt_count_UCSC'
                'i_t_ref_count_UCSC'
                'i_t_depth_UCSC'
                'i_n_depth_UCSC'
                'i_n_alt_countUCSC'
                'i_PoN_low_cov_N2913'
                'i_PoN_ok_N2913'
                'i_PoN_alt01_N2913'
                'i_PoN_alt03_N2913'
                'i_PoN_alt1_N2913'
                'i_PoN_alt3_N2913'
                'i_PoN_alt20_N2913'
                'i_PoN_germline_N2913'
                'i_PositionSampleIndex'
                'i_linearPositionSampleIndex'
                'i_PoN_low_cov_N4513'
                'i_PoN_ok_N4513'
                'i_PoN_alt01_N4513'
                'i_PoN_alt03_N4513'
                'i_PoN_alt1_N4513'
                'i_PoN_alt3_N4513'
                'i_PoN_germ3_N4513'
                'i_PoN_germ10_N4513'
                'i_PoN_Germline_N4513'
                'i_PoN_Artifact_N4513'
                'i_Variant_Classification_splice3bp'
                'i_PoN_Germline_N2913'
                'i_PoN_Artifact_N2913'
                'i_num_Centers'
                'validation_judgement_other'
                'validation_power_other'
                'validation_tumor_alt_count_other'
                'validation_tumor_ref_count_other'
                'validation_normal_alt_count_other'
                'validation_normal_ref_count_other'
                'discovery_tumor_alt_count_other'
                'discovery_tumor_ref_count_other'
                'discovery_normal_alt_count_other'
                'discovery_normal_ref_count_other'
                'validation_judgement_rna'
                'validation_power_rna'
                'validation_tumor_alt_count_rna'
                'validation_tumor_ref_count_rna'
                'validation_normal_alt_count_rna'
                'validation_normal_ref_count_rna'
                'Start_Position'
                'End_Position'
                'pon_total_counts_coverage_less_8'
                'pon_total_counts_coverage_greater_8'
                'pon_alt_count_greater1_af_greater_01percent'
                'pon_alt_count_greater2_af_greater_03percent'
                'pon_alt_count_greater3_af_greater_1percent'
                'pon_alt_count_greater3_af_greater_3percent'
                'pon_alt_count_greater3_af_greater_20percent'
                'pon_alt_count_greater10_af_greater_20percent'
                'pon_loglike'
                'pon_weight'
                'pon_pass_grid'
                'pon_pass_loglike'
                'pon_pass_weight'
                'pon_coding'
                'pon_low_alt_count'
                'i_picard_oxoQ'
                'min_val_count_rna'
                'validation_judgement_wgs'
                'validation_power_wgs'
                'validation_tumor_alt_count_wgs'
                'validation_tumor_ref_count_wgs'
                'validation_normal_alt_count_wgs'
                'validation_normal_ref_count_wgs'
                'discovery_tumor_alt_count_wgs'
                'discovery_tumor_ref_count_wgs'
                'discovery_normal_alt_count_wgs'
                'discovery_normal_ref_count_wgs'
                'min_val_count_wgs'
                'igv_bad'
                'ucsc'
                'hgsc'
                'broad'
                'bcgsc'
                'all_sites'
                'Added_by_deep_seq'
                'contEst'
                'p_ContEst'
                'ContEstKeep'
                'SCNA_NA'
                'SCNA_NB'
                'SCNA_q_hat'
                'SCNA_CCF_hat'
                'SCNA_CCF_Low'
                'SCNA_CCF_High'
                'SCNA_tau'
                'IS_SCNA'
                'purity'
                'ploidy'
                'dna_fraction_in_tumor'
                'CCF_hat'
                'CCF_Low'
                'CCF_High'
                'CCF_CI95_low'
                'CCF_CI95_high'
                'CCF_mode'
                'CCF_mean'
                'CCF_median'
                'clonal'
                'CCF_CI_low'
                'CCF_CI_high'
                
                };
            k=find(ismember(r,ft));
            cut='1';
            for k1=k(2:end)
                k1;
                cut=[cut ',' num2str(k1)];
            end
            cmd=['cut -f' cut ' ' maf_file ' > tmp.maf'];
            [o r]=unix(cmd);
            maf_file='tmp.maf'   ;
        end
    end
end
cmd =['grep "^#" ' maf_file ' | wc -l'];
[o nh]=unix(cmd);
nh=round(abs(str2num(nh)));
if ~(nh>=0)
    o
    nh
    maf_file
    keyboard
end
if nh>0
    X0 = readtable(maf_file,'Delimiter','\t','FileType','text','HeaderLines',nh);
else
    X0 = readtable(maf_file,'Delimiter','\t','FileType','text');
end    
X = table2struct(X0,'ToScalar',true);

if ~iscell(X.Chromosome)
    X.Chromosome=num2chrom(X.Chromosome)
end

function X = old_load_maf(maf_file,OPT)
% X = load_maf(maf_file,PAR)
if nargin<2
    OPT.MAF=true;
    OPT.cellstrings=true;    
end

options='';
if isfield(OPT,'cellstrings')
    options='cellstrings';
end

X=load_table(maf_file,char(9),[], [], [],options)

fn={'t_alt_count','t_ref_count','i_tumor_f','i_t_lod_fstar',...
    'i_init_n_lod','i_n_alt_count','i_n_ref_count',...
    'i_COSMIC_n_overlapping_mutations','i_ESP6500_AA_AF','i_ESP6500_EA_AF',...
    'i_1000Gp1_AF','i_1000Gp1_AFR_AF','i_1000Gp1_AMR_AF','i_1000Gp1_ASN_AF','i_1000Gp1_EUR_AF'};
fpon={'PoN_low_cov','PoN_ok','PoN_alt01','PoN_alt03','PoN_alt1','PoN_alt3','PoN_alt20','PoN_germline'};
fn=[fn fpon];
for i=1:length(fn)
    if isfield(X,fn{i})
        display([' convert ' fn{i}]);
        X.(fn{i})=str2double(X.(fn{i}));
    end
end

X.t_depth=X.t_alt_count+X.t_ref_count;
X.i_tumor_f=X.t_alt_count./X.t_depth;

% normalize PON
if all(isfield(X,fpon))
    NPON=max(X.PoN_low_cov+X.PoN_ok+X.PoN_alt01+X.PoN_alt03+X.PoN_alt1+X.PoN_alt3+X.PoN_alt20+X.PoN_germline);
    for i=1:length(fpon)
        display([' normalize ' fpon{i}]);
        X.(fpon{i})=X.(fpon{i})/NPON;
    end
end

function test

A='/Users/stewart/Projects/Cancer/Pediatric/EwingsSarcoma/PoN/filter/'
W='AN_SIGMA_EwingsSarcoma_tumorNormal'; S='AN_Sigma_Ewings_Sarcoma_07Jul2013_26TN'
S='AN_Sigma_Ewings_Sarcoma_6May2014_92T'
f1=[A W '.' S '.PoN.all.maf']
X=load_maf(f1)

function test2
X=load_maf('/Users/stewart/Projects/Cancer/DLBCL/maf/Lohr_and_Ricover_Pairs.snp.maf.annotated','trim')

