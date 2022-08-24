# Before you begin!

First, open the `settings.sh` file, and edit the very first (non-commented) line

```bash
export PROJECT_DIR="<PATH_TO_THIS_REPOSITORY>"
```

For example:

```bash
export PROJECT_DIR="/home/teg/MDSINE2_figures"
```


# Navigation

Each subdirectory is a collection of scripts for running different parts of analysis found in the paper.

All child scripts should be run from here. 
For example:
```bash
bash runtime_benchmark/initialize.sh
```

For details on running a particular piece, use the following sub-documentation:

1. [Preprocessing](preprocess/README.md)
2. [Analysis](analysis/README.md)
3. [Synthetic Benchmark](synthetic/README.md)
4. [Runtime Benchmark](runtime_benchmark/README.md)
