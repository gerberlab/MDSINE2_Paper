# Synthetic data benchmarking


# Step 1: Create the dataset

```bash
bash synthetic/create_datasets.sh
```


# Step 2: Run inference using the different methods

```bash
bash synthetic/mdsine2_infer_all.sh
bash synthetic/mdsine1_infer_all.sh
bash synthetic/run_regression_all.sh
```


# Step 3: Compute the errors, and plot them.

```bash
bash synthetic/evaluate.sh
bash synthetic/draw_plots.sh
```
