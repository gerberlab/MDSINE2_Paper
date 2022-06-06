# MDSINE2_figures


All scripts should be run from the ```scripts``` directory (to match relative pathing of 
the root `settings.sh` environment file.)

##1. Initialize the files 
-Generates the files containing information about the abundance and perturbation in a format that is compatible with code 
-The resulting files are saved in  "scripts/runtime_benchmark/other_methods/code_repo/data/gibson"
```
bash script/runtime_benchmark/initialize.sh
```

##2. Run the code to estimate the parameters of the model 
-For gLV and gLV we use both elastic-net and ridge regression; for all other models we use elastic-net 
-The outputs are saved in ../datasets/gibson/outputs (relative to ```scripts``` directory)

```
bash script/runtime_benchmark/other_methods/run_clv.sh
```
```
bash script/runtime_benchmark/other_methods/run_glv.sh
```
```
bash script/runtime_benchmark/other_methods/run_glv_ra.sh
```
```
bash script/runtime_benchmark/other_methods/run_lra.sh
```
