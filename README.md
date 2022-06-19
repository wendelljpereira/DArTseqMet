[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Computational protocol for the analysis of DArTseqMet data.

This repository contains the source code necessary to analyze DArTseqMet data, identifying the DNA methylations present in a sample on a genome-wide scale.

The method shown here was first described in the paper ["A cost-effective approach to DNA methylation detection by Methyl Sensitive DArT sequencing."](https://doi.org/10.6084/m9.figshare.10305431). Later, the same approach was used to investigate DNA methylation in clones of _Eucalyptus grandis_ planted in contrasting environments, as described in the paper ["Patterns of DNA methylation changes in elite Eucalyptus clones across contrasting environments"](https://doi.org/10.1016/j.foreco.2020.118319).

## Installation

This computational protocol is designed to be executed using the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/).

The recommended method for [installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) is by using Conda/Mamba, as shown below:

```sh
conda install -n base -c conda-forge mamba
```

Conda/mamba allows you to create different environments containing files, packages, and their dependencies that will not interact with other environments. Therefore, creating a new environment to contain the dependencies to execute this workflow is advantageous. For more information about conda and conda environments, please visit: https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

Here, we create an environment named DArTseqMet while also installing Snakemake within it. Since the workflow uses bowtie2, which depends on python2.7, we need to create a second environment that can be used by Snakemake to avoid conflits between software that depends on different versions of Python. 

```sh
conda activate base
mamba create -c conda-forge -c bioconda -n DArTseqMet snakemake

# Here we create a new environment for bowtie two, then export it as a yaml file that can be called in the rules of the worflow that depend on it.
mamba create -c conda-forge -c bioconda -n bowtie2_env python=2.7 bowtie2

conda activate bowtie2_env
conda env export > bowtie2.yaml
conda deactivate
```

Next, we activate the DArTseqMet and install other necessary software for the workflow.

```sh
conda activate DArTseqMet

## Installing R
conda install r=3.5.1
## Installing Trimmomatic from the bioconda channel
conda install -c bioconda trimmomatic
## Installing samtools from the bioconda channel
conda install -c bioconda samtools
## Installing bedtools from the bioconda channel
conda install -c bioconda bedtools
## Installing subread from the bioconda channel
conda install -c bioconda subread
## Installing fastqc from the bioconda channel
conda install -c bioconda fastqc
```

We also need to install some R packages.

```sh
conda install -c bioconda bioconductor-biostrings
conda install -c r r-docopt
```

## Executing the analysis

### Adjusting the config.yaml file.

Users need to adjust the config.yaml file to inform samples names and other parameters...

### Executing the workflow

After adjusting the config.yaml file, the execution of the worflow is as simple as running one command line.

```sh
snakemake -p -c 7 --use-conda all
```

Note that some parameters are required:

* "-c" control the number of cores to be used. This value must be the same as informed in the configuration file (config.yaml)
* "--use-conda" allows the workflow to take advantage of conda to build the environment for bowtie 2.

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.