#!/bin/bash

echo "Running all models. This will take a while to complete. In case you want to run a specific model, comment out the ones you don't want to run"

#bash cv_comparator_methods/run_lra.sh "elastic-net"
#bash cv_comparator_methods/run_glv.sh "ridge"
#bash cv_comparator_methods/run_glv_ra.sh "ridge"
#bash cv_comparator_methods/run_clv.sh "elastic-net"
bash cv_comparator_methods/run_glv.sh "elastic-net"
#bash cv_comparator_methods/run_glv-ra.sh "elastic-net"