# MDSINE2 analysis of mouse dataset


# Step 1: Fit a negative-binomial distribution for reads.

```bash
bash analysis/mdsine2_negbin.sh
```

This also outputs plots to the target directory 
`<REPO_DIR>\datasets\gibson\output\mdsine2\negbin\replicates\posterior`.


# Step 2: Run MDSINE2 core inference.

```bash
bash analysis/mdsine2_infer.sh
```


# Step 3: Run MDSINE2 inference in fixed-cluster mode.

```bash
bash analysis/mdsine2_infer_fixedcluster.sh
```

# Step 4: TODO
