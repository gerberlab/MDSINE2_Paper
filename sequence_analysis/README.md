## Making an environment for R and connecting to Jupyter

First we need to install the appropriate packages into a new `conda` environment called __rmicrobiome__.

```bash
conda create -n rmicrobiome \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  --override-channels \
  bioconductor-dada2=1.10.0

```
Now activate the new `conda` environment
```bash
conda activate rmicrobiome # this will be "source activate" on a mac
```

Install each of these one line at a time
```bash
conda install -c bioconda bioconductor-phyloseq
conda install -c bioconda r-r.utils
conda install -c r r-irkernel
conda install -c conda-forge jupyterlab
```
If for whatever reason phyloseq does not install, then go to the bottom of this markdown file


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
>[1] ‘1.10.0’

```R
packageVersion("phyloseq")
```
>[1] ‘1.26.0’

```R
quit()
```

You are now back in the terminal. Double check to make sure you have the __R-microbiome__ kernel. Refresh the jupyterhub wesbite and click on the __new__ button


![](Figures/r.png)

[//]: <> (pandoc -f markdown -t html5 -o making_r_env.html making_r_env.md --css=../css/github.css --self-contained --highlight-style=haddock --metadata pagetitle="readme")

