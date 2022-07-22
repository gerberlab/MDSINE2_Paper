# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export PROJECT_DIR="/data/cctm/darpa_perturbation_mouse_study/sawal_test/mdsine2_methods_final/MDSINE2_figures"

# ========= Don't modify below
cd ${PROJECT_DIR}/scripts
source ${PROJECT_DIR}/scripts/cv_comparator_methods/settings.sh
cd -

export LSF_GROUP=LSF_JOB_GROUP
export LSF_DIR="./lsf_files"
export LSF_QUEUE="gpu"
