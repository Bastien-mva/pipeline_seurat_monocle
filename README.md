# Cell specialization and coordination in Arabidopsis leaves upon pathogenic attack revealed by scRNA-seq

The code for the paper [Cell specialization and coordination in Arabidopsis
leaves upon pathogenic attack revealed by scRNA-seq](https://www.biorxiv.org/content/10.1101/2023.03.02.530814v1). The data is available
[here](link). You should put the downloaded data inside this respository to run
the script ```pipeline.R```. Many of the available plots in the paper are taken
from this script, and all the plots coming from Seurat and Monocle can be
derived changing only arguments in the ```pipeline.R``` (for example changing
the name of genes for UMAP plots).

To run the script, you must install the following packages:

- dplyr
- Seurat
- ggplot2
- patchwork
- ggpubr
- monocle3
- SeuratWrappers

Once installed, you can run ```r pipeline.R ```

## get_citation


