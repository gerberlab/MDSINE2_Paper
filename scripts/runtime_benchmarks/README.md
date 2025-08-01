NOTE: This is NOT the analysis that generated Supp. Figure 7 (runtime analysis). That was done using a newer script (`semisynthetic2/inference/runtime_benchmark.sh`)



Run the scripts in the following order:

1. `runtime_benchmarks/mdsine2_filter.sh`

2. `runtime_benchmarks/mdsine2_run_replicates.sh`

3. `runtime_benchmarks/mdsine2_infer.sh`

Each script depends on the previous.
Script #1 depends on `preprocess` pipeline having finished. In particular, it looks for `gibson_healthy_agg_filtered.pkl`
in the preprocessed directory.