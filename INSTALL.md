# Snapshot Installation Guide

## Dependencies
The following packages and tools are required to run Snapshot:
- Python3:
  - numpy
  - sklearn
- R:
  - ggplot2
  - pheatmap
  - igraph
  - networkD3
  - data.table
  - mclust
  - dplyr
  - lsa
  - cba
  - RColorBrewer
  - tidyverse
- bedtools
- ucsc_tools


## Installation Steps
1. Clone the Github repository:
```
git clone https://github.com/guanjue/snapshot.git
```
2. Update conda channels
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek
```
3. Set up a conda environment named "snapshot":
```
conda create -n snapshot python=3 r=3.6 bedtools ucsc_tools numpy scikit-learn r-ggplot2 r-pheatmap r-igraph r-networkD3 r-data.table r-mclust r-dplyr r-lsa r-cba r-RColorBrewer r-tidyverse
conda activate snapshot
```

Note: Detailed instructions on how to install conda can be found in the [conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html).
