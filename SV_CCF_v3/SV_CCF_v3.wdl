task sv_ccf_v3_task_1 {
    File SV_BP_FILE
    File ACS_CCF_FILE
    String GENDER
    String ID

    Float? ram_gb
    Int? local_disk_gb
    Int? num_preemptions

    command {
        set -euo pipefail
     pwd
        ls -lath /src
        ls -lath 

        echo "SV_BP_FILE: ${SV_BP_FILE}"
        echo "ACS_CCF_FILE: ${ACS_CCF_FILE}"
        echo "GENDER: ${GENDER}"
        echo "ID: ${ID}"

        /src/run_fc_SV_CCF_v3.sh /usr/local/MATLAB/MATLAB_Runtime/v901 ${SV_BP_FILE} ${ACS_CCF_FILE} ${GENDER}  ${ID} 
        ls -lath 
    }

    output {
          File SV_BP_CCF_v3_tsv="${ID}.SV_CCF_v3.tsv"
    }

    runtime {
        docker : "chipstewart/sv_ccf_v3_task_1:1"
        memory: "${if defined(ram_gb) then ram_gb else '2'}GB"
        disks : "local-disk ${if defined(local_disk_gb) then local_disk_gb else '10'} HDD"
        preemptible : "${if defined(num_preemptions) then num_preemptions else '0'}"
    }

    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }
}

workflow SV_CCF_v3 {
    call sv_ccf_v3_task_1 
}
