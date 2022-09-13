# submits the jobs for running cLV model
set -e
source settings.sh


bash create_regression_job.sh "clv" "elastic-net"
bash create_regression_job.sh "glv" "ridge"
bash create_regression_job.sh "glv" "elastic-net"
bash create_regression_job.sh "glv-ra" "ridge"
bash create_regression_job.sh "glv-ra" "elastic-net"
bash create_regression_job.sh "lra" "elastic-net"
bash create_mdsine2_job.sh 2
bash create_mdsine2_job.sh 3
bash create_mdsine2_job.sh 4
bash create_mdsine2_job.sh 5
bash create_mdsine2_nomodules_job.sh 2
bash create_mdsine2_nomodules_job.sh 3
bash create_mdsine2_nomodules_job.sh 4
bash create_mdsine2_nomodules_job.sh 5


run_script_path=${LSF_DIR}/submit_all.sh

cat <<- EOFDOC > $run_script_path
#!/bin/bash
echo "Searching for LSF jobs to run."
for f in */*.lsf; do
	echo "Submitting job: \$f"
	bsub < \$f
done
EOFDOC
