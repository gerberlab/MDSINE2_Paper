# MDSINE2_Paper

Pull this code via
```bash
git clone --recurse-submodules https://github.com/gerberlab/MDSINE2_Paper.git
```


# Important links
- Associated GitHub repo for the Bayesian model: ["MDSINE2"](https://github.com/gerberlab/MDSINE2)<br />
  <a href="https://github.com/gerberlab/MDSINE2"><img alt="" src="https://img.shields.io/badge/GitHub-MDSINE2-blue?style=flat&logo=github"/></a>
- Raw sequences from longitudinal experiments on NCBI: <br />
- <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA784519/"><img alt="" src="https://img.shields.io/badge/NCBI-PRJNA784519-blue?style=flat"/></a>
- Zenodo repositories containing pre-computed outputs of inference of the full dataset and cross-validation: <br />
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8006853.svg)](https://doi.org/10.5281/zenodo.8006853) (CV + Replicates) <br />
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8208503.svg)](https://doi.org/10.5281/zenodo.8208503) (Full Inference) <br />
  (Note: please ensure downloads are from the latest versions of these records. To extract a `.tar.xf` archive on Linux, use the `tar xf` command; on Windows, use 7-zip.)


# Jupyter Notebooks
- Folder containing [tutorials as notebooks exploring the model, data and paper](https://github.com/gerberlab/MDSINE2_Paper/tree/master/tutorials) that can be opened directly in Google Colab:<br />
<a href="https://github.com/gerberlab/MDSINE2_Paper/tree/main/tutorials"><img alt="" src="https://img.shields.io/badge/Jupyter Notebooks-MDSINE2%20Tutorials-blue?style=flat&logo=jupyter"/></a>
- Folder containing [notebooks to reproduce figures found in the paper](https://github.com/gerberlab/MDSINE2_Paper/tree/master/paper_figures) (requires local Jupyter instance to re-run)


# Instructions

For details about individual steps in the analysis, refer to the primary internal README: [scripts/README.md](scripts/README.md)

As noted in that document, all scripts should be run from the `scripts` directory.
This requirement is in place to match relative pathing of the root `scripts/settings.sh` environment file, from which all other paths are defined relatively.

## Setup for analysis on a local machine

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

Once the above installation is done, one should also install Jupyter/Jupyterlab for data visualization.
```
conda install -c conda-forge jupyter jupyterlab ipywidgets nodejs
jupyter nbextension enable --py widgetsnbextension
```

## Analysis on a local machine

We provided shell scripts that implement several pipelines. Link: [Analysis Pipeline](scripts/README.md) 
 
To reproduce the figures using the output of these analyses, please go through the Jupyter notebooks 
located in the `paper_figures` directory of this repository.
To run them yourself, you will need to run a jupyter instance, for example:
```bash
cd MDSINE2_Paper
jupyter lab --port 8888
```
and navigate to `paper_figures/` through JupyterLab.
