Order of scripts:

1) pick_ground_truth.sh - evaluates real-data inference pickle files and picks a "best" GLV ground truth paramset.
2) create_initial_conditions.sh - generates all initial conditions for all trajectory replicates. 
3) generate_trajectories.sh - generates all fwsim trajectories for all trajectory replicates, for all perturbation counts.
4) [branching points:]
   a) create_dataset_tables_vary_mice.sh - 
   b) create_dataset_tables_vary_perts.sh - 
   c) create_dataset_tables_thin_timepoints.sh - 