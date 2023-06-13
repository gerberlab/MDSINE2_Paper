# MDSINE2 - Full analysis of mouse dataset

Prerequisites: These scripts assume that the full pipeline found in `preprocess` has been run.

# Step 1: Fit a negative-binomial distribution for reads.

```bash
bash analysis/mdsine2_negbin.sh
```

This outputs inference pickle files (and accompanying plots) into
`MDSINE2_Paper\datasets\gibson\output\mdsine2\negbin\replicates`.


# Step 2: Run MDSINE2 core inference

Next, run inference using ten multiple seeds. 
Inference is programmed by this script to run in series (one after the other).

```bash
bash analysis/mdsine2_infer_all.sh
```

This script will automatically:
1. Invoke `mdsine2 infer` on ten different seeds (using the settings found in `settings.sh`),
2. Plot each seed's posterior using `mdsine2 visualize-posterior`, and
3. Merge all seeds' posteriors, by dumping samples into concatenated numpy arrays using `mdsine2 extract-posterior`.

All outputs are located in `MDSINE2_Paper/datasets/gibson/output/mdsine2/inference`.
Each seed is located in the respective subdirectory `healthy-seed<SEED>`,
and the merged files are dumped into separate numpy arrays (for each gLV parameter) 
in the `merged_studies` subdirectory.

# Step 3: Run MDSINE2 inference in fixed-cluster mode.

Next, run fixed-cluster inference using the concatenated seeds' posterior-approximating markov chains.

```bash
bash analysis/mdsine2_infer_fixedcluster.sh
```

Outputs are located in `MDSINE2_Paper/datasets/gibson/output/mdsine2/inference/merged_studies_fixed_cluster`.

# Optional: evaluate keystoneness metrics.

This is required only if rendering the keystoneness plot (Figure 5).

```bash
bash analysis/keystoneness.sh
```

Outputs are located in `MDSINE2_Paper/datasets/gibson/output/mdsine2/inference/merged_studies/keystoneness`.
