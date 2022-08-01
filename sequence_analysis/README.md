Install an R environement

conda install -n envname -c r r
conda activate envname



R

install.packages("BiocManager")
BiocManager::install("dada2")
