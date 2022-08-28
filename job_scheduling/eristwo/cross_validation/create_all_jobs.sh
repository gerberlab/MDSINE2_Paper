# submits the jobs for running cLV model
set -e
source settings.sh


bash create_regression_job.sh "clv" "elastic-net"
bash create_regression_job.sh "glv" "ridge"
bash create_regression_job.sh "glv" "elastic-net"
bash create_regression_job.sh "glv-ra" "ridge"
bash create_regression_job.sh "glv-ra" "elastic-net"
bash create_regression_job.sh "lra" "elastic-net"
