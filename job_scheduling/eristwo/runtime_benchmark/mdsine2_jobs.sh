#!/bin/bash

# RUN THIS SCRIPT FROM job_scheduling/eristwo/runtime_benchmark/
set -e
source settings.sh


lsf_subdir=${LSF_DIR}/mdsine2
log_dir=${lsf_subdir}/logs
mkdir -p ${lsf_subdir}
mkdir -p ${log_dir}


generate_lsf()
{
	n_taxa=$1
	trial=$2

	inference_lsf_path="${lsf_subdir}/taxa_${n_taxa}_trial_${trial}.lsf"

	echo "Creating ${inference_lsf_path}"
	cat <<- EOFDOC > $inference_lsf_path
#!/bin/bash
#BSUB -J taxa_${n_taxa}_trial_${trial}
#BSUB -g ${LSF_JOB_GROUP}/mdsine2
#BSUB -o ${log_dir}/filter_reads_${n_reads}_qs_${quality_shift}_trial_${trial}-%J.out
#BSUB -e ${log_dir}/filter_reads_${n_reads}_qs_${quality_shift}_trial_${trial}-%J.err
#BSUB -q ${LSF_QUEUE}
#BSUB -n ${LSF_N_CORES}
#BSUB -M ${LSF_MEM}
#BSUB -R rusage[mem=${LSF_MEM}]

cd ${PROJECT_DIR}/scripts
bash runtime_benchmark/mdsine2/mdsine2_infer.sh $n_taxa $trial
EOFDOC
}


for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${MDSINE2_NUM_TRIALS}+1; trial++ )); do
		generate_lsf $n_taxa $trial
	done
done
