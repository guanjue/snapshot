conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek

conda create -n snapshot python=3 r=3.6 bedtools ucsc_tools numpy scikit-learn r-ggplot2 r-pheatmap r-igraph r-networkD3 r-data.table r-mclust r-dplyr r-lsa r-cba r-RColorBrewer r-tidyverse
conda activate snapshot
