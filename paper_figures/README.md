# To run these notebooks

Download all files from Zenodo, or run inference from scratch.
If extracting from Zenodo, files should be placed in the same manner as inference would be done from scratch. 
This means that the following directory structure (and associated files inside the directories) must be present:

1) Preprocessed inputs:
    - `MDSINE2_Paper/datasets/gibson/healthy/preprocessed` -- must contain the following input files:
        - gibson_healthy_agg_filtered.pkl
        - gibson_replicates_agg_filtered.pkl
        - gibson_inoculum_agg_filtered.pkl
3) Main inference:
    - `MDSINE2_Paper/datasets/gibson/healthy/output/mdsine2` -- must contain the following folders which contain outputs from MDSINE2.
        - inference/healthy-seed[10, 11, 12, 13, 14, 15, 16, 17, 19, 20]
        - inference/merged_studies
        - inference/merged_studies_fixed_cluster
    - `MDSINE2_Paper/datasets/semisynthetic2` -- must contain the following subdirs:
        - trajectory_replicate_[0-9]
        - truth
4) Cross-validation inference:
    - `MDSINE2_Paper/datasets/gibson/healthy/cross_validation` -- must contain the following folders, which contain the cross-validation outputs.
        - mdsine2-modules
        - mdsine2-nomodules
        - regression_clv_elastic_net
        - regression_glv_elastic_net
        - regression_glv-ra_elastic_net
        - regression_glv-ra_ridge
        - regression_glv_ridge
        - regression_lra_elastic_net
        - edge_density_[1-5]
        - growth_si_var_[1-4]
        - interaction_var_[1-4]
        - pert_var_[1-4]
    - `MDSINE2_Paper/datasets/gibson/uc/cross_validation` -- must contain the following folders, which contain the cross-validation outputs.
        - mdsine2-modules
        - mdsine2-nomodules
        - regression_clv_elastic_net
        - regression_glv_elastic_net
        - regression_glv-ra_elastic_net
        - regression_glv-ra_ridge
        - regression_glv_ridge
        - regression_lra_elastic_net
     
# TODO
1) edit fig4_semisynthetic_v2_cache to point to local dir, instead of cctm briefcase.
2) edit sup_cv_prior_comparison to point to local dir, instead of cctm briefcase.
3) test sup_cv_prior_comparison to load from default run, instead of default_new.

# Notebook file list

Notebooks for final paper submission:

fig2_mice.ipynb -- Main figure 2 (Mice experiment and taxa visualization)
fig3_cross_validation.ipynb -- Main figure 3 (Cross validation experiment, render healthy and uc cohorts)
    -> Latest version is fig4_cross_validation-BOTH.ipynb (both healthy/uc cohorts in one figure)
    -> Renamed to fig4_cross_validation.ipynb (Jul 2, 2025)
    -> Renamed to fig3_cross_validation.ipynb (July 14, 2025)
fig4_semisynthetic_v2_cache.ipynb -- Semisynthetic benchmark (Main figure 4)
fig5_modules.ipynb -- OTU breakdown by modules + abundance heatmap.
fig6_stability.ipynb -- Eigenvalue calculations of interaction matrices.

sfig_small_synthetic.ipynb -- Small synthetic experiment (Supp. Fig 10)
sfig_enrichment.ipynb -- per-module taxonomic composition/enrichment. (Supp. Fig 4)