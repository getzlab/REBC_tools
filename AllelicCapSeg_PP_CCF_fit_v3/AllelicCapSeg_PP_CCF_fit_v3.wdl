task alleliccapseg_pp_ccf_fit_v3_task_1 {
    
    File ACS_FILE
    String PURITY
    String PLOIDY
    String ID

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions
    

    command {
        set -euo pipefail
        pwd
        ls -lath /src
        ls -lath 

        echo "FILE: ${ACS_FILE}"
        echo "PURITY: ${PURITY}"
        echo "PLOIDY: ${PLOIDY}"
        echo "ID: ${ID}"

        /src/run_fc_AllelicCapSeg_PP_CCF_fit_v3.sh /usr/local/MATLAB/MATLAB_Runtime/v901 ${ACS_FILE} ${PURITY} ${PLOIDY} ${ID} female
        ls -lath 
    }

    output {
        File alleliccapseg_pp_ccf_fit_v3_out="${ID}.acs.ccf.tsv"
        File alleliccapseg_pp_ccf_fit_v3_collapsed_out="${ID}.acs.ccf.collapse.tsv"
    }

    runtime {
        docker : "chipstewart/alleliccapseg_pp_ccf_fit_v3_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow AllelicCapSeg_PP_CCF_fit_v3 {

    call alleliccapseg_pp_ccf_fit_v3_task_1 
}
