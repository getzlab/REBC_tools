#Remove Firecloud filtered calls with less than three callers and more than zero alt reads in normal:
  
SVNALG=/GoogleDrive/Cancer/REBC/SV/REBC_primary_pair_393.SV_BP_CCF_v3.noHot.04May2020.tsv
SVFILT=/GoogleDrive/Cancer/REBC/SV/REBC_primary_pair_393.SV_BP_CCF_v3.noHot.24May2020.FILT.tsv
python dedup_svs_v3.py -i $SVNALG -s 350 -N 0 -A 3 -v 0.05 -T 3 > $SVFILT
wc -l $SVFILT

#Produces a filtered SV list with 1139 calls, removing 35 SV with VCF_NALT>0 and NALG<3

#Pre-filtered Fireclould calls (no NALG, TALT, VAF, blacklist, or hotspot filters) 
SVNALG=/GoogleDrive/Cancer/REBC/SV/REBC_primary_pair_393.aggregated.all.NALG.tsv
SVFILT=/GoogleDrive/Cancer/REBC/SV/REBC_primary_pair_393.aggregated.all.filt.24May2020.tsv
python scripts/dedup_svs_v3.py -i $SVNALG -s 350 -N 2 -A 3 -v 0.075 -T 2 > $SVFILT
wc -l $SVFILT

#Produces a SV list with 9810 calls, removing 76740 SV with (VCF_NALT=1 and NALG<3) or (VAF<-0.075) or (VCF_TALT<2)
