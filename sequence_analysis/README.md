# Processing raw reads from sequencing and turning them into ASVs


## Making an environment for R and connecting to Jupyter (for convenience if you want to turn the R script into a notebook later)

First we need to install the appropriate packages into a new `conda` environment called __rmicrobiome__.

```bash
conda create -n rmicrobiome \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  --override-channels \
  bioconductor-dada2=1.20

```
Now activate the new `conda` environment
```bash
conda activate rmicrobiome
```

Install each of these one line at a time
```bash
conda install -c bioconda bioconductor-phyloseq
conda install -c bioconda r-r.utils
conda install -c r r-irkernel
conda install -c conda-forge jupyterlab
```

Now we are linking this environment to `Jupyter`.

```bash
R -e 'IRkernel::installspec(name = "rmicrobiome", displayname = "R-microbiome")'
```
>[InstallKernelSpec] Installed kernelspec rmicrobiome in /jupyteruser/50940993/.local/share/jupyter/kernels/rmicrobiome

Open `R`
```bash
R
```
Check package versions to ensure things were installed properly

```R
packageVersion("dada2")
```
>[1] ‘1.20.0’

```R
packageVersion("phyloseq")
```
>[1] ‘1.26.0’


## Other files you will need in place

Now you are ready to run the main script `dada2.R`
