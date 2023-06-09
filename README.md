# MDSINE2_figures

Pull this code via
```bash
git clone --recurse-submodules https://github.com/gibsonlab/MDSINE2_figures.git
```

**Abstract**: Despite the importance of microbial dysbiosis in human disease, the phenomenon remains poorly understood. We provide the first comprehensive and predictive model of dysbiosis at ecosystem-scale, leveraging our new machine learning method for efficiently inferring compact and interpretable dynamical systems models. Coupling this approach with the most densely temporally sampled interventional study of the microbiome to date, using microbiota from healthy and dysbiotic human donors that we transplanted into mice subjected to antibiotic and dietary interventions, we demonstrate superior predictive performance of our method over state-of-the-art techniques. Moreover, we demonstrate that our approach uncovers intrinsic dynamical properties of dysbiosis driven by destabilizing competitive cycles, in contrast to stabilizing interaction chains in the healthy microbiome, which have implications for restoration of the microbiome to treat disease.

Important links
- Associated GitHub repo for the Bayesian model: ["MDSINE2"](https://github.com/gerberlab/MDSINE2)<br />
  <a href="https://github.com/gerberlab/MDSINE2"><img alt="" src="https://img.shields.io/badge/GitHub-MDSINE2-blue?style=flat&logo=github"/></a>
- Folder containing [tutorials as notebooks exploring the model, data and paper](https://github.com/gerberlab/MDSINE2_Paper/tree/master/google_colab) that can be opened directly in Google Colab<br />
<a href="https://github.com/gerberlab/MDSINE2_Paper/tree/main/tutorials"><img alt="" src="https://img.shields.io/badge/Jupyter Notebooks-MDSINE2%20Tutorials-blue?style=flat&logo=jupyter"/></a>
- Raw sequences from longitudinal experiments on NCBI <br />
<a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA784519/"><img alt="" src="https://img.shields.io/badge/NCBI-PRJNA784519-blue?style=flat"/></a>


# Instructions

For details about individual steps in the analysis, refer to the primary internal README: [scripts/README.md](scripts/README.md)

As noted in that document, all scripts should be run from the `scripts` directory.
This requirement is in place to match relative pathing of the root `scripts/settings.sh` environment file, from which all other paths are defined relatively.

## Setup: Analysis on a local machine

This section outlines how to run MDSINE2 analysis on our dataset in full, with `bash`, `conda` and `git`.

One must first install the MDSINE2 package.
The recommended setup uses the conda recipe provided in that repository.

```
cd ~/work
git clone https://github.com/gerberlab/MDSINE2
cd MDSINE2
conda env create -f conda_recipe.yml 
conda activate mdsine2
```

Next, clone this repository which contains the data and scripts to perform the analysis.

```
cd ~/work
git clone https://github.com/gerberlab/MDSINE2_Paper
cd MDSINE2_Paper
```

Once the above installation done, one should also install some additional modules, including Jupyter.
```
conda install -c conda-forge jupyterlab
```
