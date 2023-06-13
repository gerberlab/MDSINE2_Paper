# Synthetic data benchmarking


# Step 1: Create the dataset

```bash
bash synthetic/small/create_datasets.sh
```


# Step 2: Run inference using the different methods

```bash
bash synthetic/small/mdsine2_infer_all.sh
bash synthetic/small/mdsine1_infer_all.sh
bash synthetic/small/run_regression_all.sh
```


# Step 3: Compute the metrics. (For plotting the data, refer to the accompanying jupyter notebook.)

```bash
bash synthetic/small/evaluate.sh
```
