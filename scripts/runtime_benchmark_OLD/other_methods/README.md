# MDSINE2_figures


All scripts should be run from the ```scripts``` directory (to match relative pathing of 
the root `settings.sh` environment file.)

## 1. Initialize the files 

-Generates the files containing information about the abundance and perturbation in a format that is compatible with code.

-The location of the source files are provided in the "settings.sh". Modify the variables to provide the right input data.

```
bash runtime_benchmark/initialize.sh
```

## 2. Run the code to estimate the parameters of the model 

-For gLV and gLV we use both elastic-net and ridge regression; for all other models we use elastic-net

```
bash runtime_benchmark/other_methods/run_clv.sh
```
```
bash runtime_benchmark/other_methods/run_glv.sh
```
```
bash runtime_benchmark/other_methods/run_glv_ra.sh
```
```
bash runtime_benchmark/other_methods/run_lra.sh
```
