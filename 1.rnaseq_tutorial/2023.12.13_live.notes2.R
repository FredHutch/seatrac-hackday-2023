library(tidyverse)
library(Seurat)

load("data/dat_singlecell.RData")
str(dat_sc)

# UMAP
# DO NOT RUN
dat_sc <- RunUMAP(dat_sc, dims = 1:10)

DimPlot(dat_sc, reduction="umap",
        group.by = "orig.ident")
DimPlot(dat_sc, reduction="umap",
        group.by = "orig.clusters")

umap_xy <- as.data.frame(dat_sc[["umap"]]@cell.embeddings)
umap_xy %>% head

cell_clusters <- as.data.frame(dat_sc$orig.clusters)
cell_clusters %>% head

umap_xy_clusters <- umap_xy %>% 
  #move rownames to data column
  rownames_to_column("cell") %>% 
  #merge with cluster info
  full_join(
    cell_clusters %>% rownames_to_column("cell")
  ) %>% 
  #clean column names
  rename(ident = `dat_sc$orig.clusters`)
umap_xy_clusters %>% head

#highlight genes
FeaturePlot(dat_sc, reduction="umap",
            features =c('STAT1','TNFRSF8'))

#tSNE
#DO NOT RUN
dat_sc <- RunTSNE(dat_sc, dims = 1:10)

DimPlot(dat_sc, reduction="tsne",
        group.by = "orig.ident")
DimPlot(dat_sc, reduction="tsne",
        group.by = "orig.clusters")
