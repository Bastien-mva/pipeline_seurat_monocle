# Cell specialization and coordination in Arabidopsis leaves upon pathogenic attack revealed by scRNA-seq


The code for the paper [Cell specialization and coordination in Arabidopsis
leaves upon pathogenic attack revealed by scRNA-seq](https://www.biorxiv.org/content/10.1101/2023.03.02.530814v1). The data is available
[here](link). First, you should clone the repo:
```
git clone https://github.com/Bastien-mva/pipeline_seurat_monocle.git
cd pipeline_seurat_monocle
```

You should put the downloaded data inside this respository. Then unzip it.
To run the script, you must install the following packages:

- dplyr
- Seurat
- ggplot2
- patchwork
- ggpubr
- monocle3
- SeuratWrappers
- ff


You  can then run
```
Rscript pipeline.R
```
Alternatively you can run it in Rstudio.

Many of the available plots in the paper are taken
from this script, and all the plots coming from Seurat and Monocle can be
derived changing only arguments in the ```pipeline.R``` script (for example changing
the name of genes for UMAP plots). Note that the script choose the root node
for pseudotime analysis automatically, so that the script can be run from the command line.
Comments indicate how to choose the root node manually if you are in
interactive mode (for example Rstudio).

## References

         Trapnell C. et. al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat. Biotechnol. 32, 381-386 (2014). https://doi.org/10.1038/nbt.2859 Qiu, X. et. al. Reversed graph embedding resolves complex single-cell trajectories. Nat. Methods 14, 979-982 (2017). https://doi.org/10.1038/nmeth.4402 Cao, J. et. al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496-502 (2019). https://doi.org/10.1038/s41586-019-0969-x
         Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reportsvolume 9, Article number: 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z Levine, J. H. et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047 Levine, J. H., et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.04


