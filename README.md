[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Computational protocol for the analysis of DArTseqMet data.

This repository contains the source code necessary to analyze DArTseqMet data, identifying the DNA methylations present in a sample on a genome-wide scale.

The method described here was first described in the paper ["A cost-effective approach to DNA methylation detection by Methyl Sensitive DArT sequencing."](https://doi.org/10.6084/m9.figshare.10305431). Later, the same approach was used to investigate DNA methylation in clones of _Eucalyptus grandis_ planted in contrasting environments, as described in the paper ["Patterns of DNA methylation changes in elite Eucalyptus clones across contrasting environments"](https://doi.org/10.1016/j.foreco.2020.118319).

## Installation

This computational protocol is designed to be executed using the [Snakemake workflow management system](https://snakemake.readthedocs.io/en/stable/).

The recommended method for [installing Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) is by using Conda/Mamba, as shown below:

```sh
conda install -n base -c conda-forge mamba
```

Conda/mamba allows you to create different environments containing files, packages, and their dependencies that will not interact with other environments. Therefore, creating a new environment to contain the dependencies to execute this workflow is advantageous. For more information about conda and conda environments, please visit: https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

Here, we create an environment named DArTseqMet while also installing Snakemake within it.

```sh
conda activate base
mamba create -c conda-forge -c bioconda -n DArTseqMet snakemake

```

```sh
conda activate DArTseqMet
```

Next, we install all dependencies within the environment we created.

```sh
## Install R
conda install r=3.5.1
## Install Trimmomatic from the bioconda channel
conda install -c bioconda trimmomatic
## Install samtools from the bioconda channel
conda install -c bioconda samtools
## Install bowtie2 from the bioconda channel
conda install -c bioconda bowtie2
## Install bedtools from the bioconda channel
conda install -c bioconda bedtools
## Install subread from the bioconda channel
conda install -c bioconda subread
## Install fastqc from the bioconda channel
conda install -c bioconda fastqc
```

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.