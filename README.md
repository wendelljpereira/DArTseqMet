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


# Bowtie2 needs python 2.7, that why we create a new environment that can be used by snakemake later
mamba create -c conda-forge -c bioconda -n bowtie2_env python=2.7 bowtie2

# exports the env as a yaml file, so snakemake can acess it in the rules where python2.7 and bowtie2 are required
conda activate bowtie2_env
conda env export > bowtie2.yaml
conda deactivate
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

We also need to install R packages.

```sh
conda install -c bioconda bioconductor-biostrings

conda install -c r r-docopt
```

## List of softwares, and their versions, used to generate the results described in the paper:

* [R](https://cran.r-project.org/) : **3.5.1**  
* [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) : **0.36**  
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) : **2.3.5**  
* [samtools](http://samtools.sourceforge.net/) : **1.8**  
* [bedtools](https://bedtools.readthedocs.io/en/latest/) : **2.27.1**  
* [subread](http://subread.sourceforge.net/) : **1.6.2**  
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) : **0.11.4**  

## List of R packages, and their versions, used to generate the results described in the paper:

* [alluvial](https://cran.r-project.org/web/packages/alluvial/index.html): **0.1.2**  
* [berryFunctions](https://cran.r-project.org/web/packages/berryFunctions/index.html): **1.18.2**  
* [corrplot](https://cran.r-project.org/web/packages/corrplot/index.html): **0.84**  
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html): **1.12.2**  
* [docopt](https://cran.r-project.org/web/packages/docopt/index.html): **0.6.1**  
* [doMC](https://cran.r-project.org/web/packages/doMC/index.html): **1.3.5**  
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html): **0.8.3**  
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html): **1.4.4**  
* [gdata](https://cran.r-project.org/web/packages/gdata/index.html): **2.18.0**  
* [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html): **0.4.6**  
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html): **3.2.1**  
* [ggsci](https://cran.r-project.org/web/packages/ggsci/index.html): **2.9**  
* [ggthemes](https://cran.r-project.org/web/packages/ggthemes/index.html): **4.2.0**  
* [googleVis](https://cran.r-project.org/web/packages/googleVis/index.html): **0.6.3**  
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html): **2.3**  
* [gridGraphics](https://cran.r-project.org/web/packages/gridGraphics/index.html): **0.4.1**  
* [Matching](https://cran.r-project.org/web/packages/Matching/index.html): **4.9.6**  
* [plyr](https://cran.r-project.org/web/packages/plyr/index.html): **1.8.4**  
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html): **1.1.2**  
* [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html): **1.4.3**  
* [sjstats](https://cran.r-project.org/web/packages/sjstats/index.html): **0.17.6**  
* [splitstackshape](https://cran.r-project.org/web/packages/splitstackshape/index.html): **1.4.8**  
* [stringr](https://cran.r-project.org/web/packages/stringr/index.html): **1.4.0**  
* [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html): **1.0.0**  
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html): **1.2.1**  
* [VennDiagram](https://cran.r-project.org/web/packages/VennDiagram/index.html): **1.6.20**  
* [AnnotationForge](https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html): **1.24.0**  
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html): **2.38.0**  
* [clusterProfiler](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html): **3.10.1**  
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): **1.22.2**  
* [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html): **3.24.3**  
* [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html): **1.34.0**  
* [ggbio](https://bioconductor.org/packages/release/bioc/html/ggbio.html): **1.30.0**  
* [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html): **1.26.5**  
* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html): **1.42.2**  



<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.