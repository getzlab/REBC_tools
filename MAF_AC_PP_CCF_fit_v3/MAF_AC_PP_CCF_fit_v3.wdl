task maf_ac_pp_ccf_fit_v3_task_1 {

    File MAF_FILE
    File ACS_FILE
    String PURITY
    String PLOIDY
    String GENDER
    String ID

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    #**Define additional inputs here**

    command {
        set -euo pipefail
    
        pwd
        ls -lath /src
        ls -lath 

        echo "MAF_FILE: ${MAF_FILE}"
        echo "ACS_FILE: ${ACS_FILE}"
        echo "PURITY: ${PURITY}"
        echo "PLOIDY: ${PLOIDY}"
        echo "GENDER: ${GENDER}"
        echo "ID: ${ID}"

        /src/run_fc_MAF_AC_PP_CCF_fit_v3.sh /usr/local/MATLAB/MATLAB_Runtime/v901 ${MAF_FILE} ${ACS_FILE} ${PURITY} ${PLOIDY} ${ID} ${GENDER}
        ls -lath 
    }

    output {
        File alleliccapseg_pp_ccf_fit_v3_out="${ID}.acs.ccf.tsv"
        File alleliccapseg_pp_ccf_fit_v3_collapsed_out="${ID}.acs.ccf.collapse.tsv"
        File alleliccapseg_pp_ccf_fit_v3_ABSOLUTE_SEGTAB="${ID}.AllelicCapSeg_PP_CCF_fit_v3.ABSOLUTE_SEGTAB.tsv"
        File maf_ac_pp_ccf_fit_v3_maf="${ID}.MAF_AC_PP_CCF_fit_v3.maf"
        File maf_ac_pp_ccf_fit_v3_ABSOLUTE_maf="${ID}.MAF_AC_PP_CCF_fit_v3.ABSOLUTE.maf"  
    }

    runtime {
        docker : "chipstewart/maf_ac_pp_ccf_fit_v3_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '7'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow MAF_AC_PP_CCF_fit_v3 {
    call maf_ac_pp_ccf_fit_v3_task_1
}
