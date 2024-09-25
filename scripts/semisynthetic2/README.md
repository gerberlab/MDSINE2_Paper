All scripts must be run from the scripts/ root directory (e.g. `bash semisynthetic2/data_generation/pick_ground_truth.sh`)

Order of scripts:

1) data_generation/pick_ground_truth.sh - evaluates real-data inference pickle files and picks a "best" GLV ground truth paramset.
2) data_generation/create_initial_conditions.sh - generates all initial conditions for all trajectory replicates. 
3) data_generation/generate_trajectories.sh - generates all fwsim trajectories for all trajectory replicates, for all perturbation counts.
4) data_generation/create_dataset_tables_raw.sh - Creates the *complete* master dataset tables. The "counts.tsv" file generated is NOT the proper format for MDSINE2; the downstream scripts will handle this. 
5) [branching points:]
   a) data_generation/create_dataset_tables_vary_mice.sh - Slice the *.tsv files from the previous step, and slices according to the desired # of mice.
   b) data_generation/create_dataset_tables_vary_perts.sh - 
   c) data_generation/create_dataset_tables_thin_timepoints.sh - 