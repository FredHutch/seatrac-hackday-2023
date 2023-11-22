library(tidyverse)
library(data.table)
library(Seurat)

# Load counts data
## This file was downloaded from GEO and not included on GitHub due to size
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139598
count <- data.table::fread("data/GSE139598_Week25_Unstimulated.Cells.UMI.collapsed.csv.gz") %>% 
  column_to_rownames("V1") %>% 
  as.matrix()

# Convert to seurat object
dat <- CreateSeuratObject(counts = count, project = "Darrah",
                           min.cells = 3, min.features = 200)
rm(count)
gc()

# Calculate % mitochondrial
## No counts present ##
# dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")

# Log normalize
dat <- NormalizeData(dat, normalization.method = "LogNormalize")

# Find highly variable features
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

#scale data
dat <- ScaleData(dat, features = rownames(dat))

#Perform PCA
dat <- RunPCA(dat, features = VariableFeatures(object = dat))

# Add original analysis clusters
clust <- read_tsv("../darrah_etal/updated_alexandria_metadata.txt") %>% 
  filter(NAME != "TYPE") %>% 
  separate(NAME, into=c("week","NAME"), sep="_", extra="merge") %>% 
  filter(week == "Week25")

clusters <- clust$cell_type__ontology_label
names(clusters) <- clust$NAME
dat[["orig.clusters"]] <- clusters

# Calculate ordinations 
dat <- RunUMAP(dat, dims = 1:10)
dat <- RunTSNE(dat, dims = 1:10)

# Down-sample to 0.5% or 10 cells per cell type
cells_keep <- c()
for(cell in unique(clust$cell_type__ontology_label)){
  temp <- clust %>% 
    filter(cell_type__ontology_label==cell)
  
  n <- round(nrow(temp)*0.005)
  if(n<10){ n <- 10 }
  
  set.seed(42)
  temp2 <- sample(temp$NAME, n)
  cells_keep <- c(cells_keep, temp2)
}

## Filter cells
dat_sub <- dat[, cells_keep]

# Full UMAP
p1 <- DimPlot(dat, reduction = "umap", group.by = "orig.clusters", 
              pt.size=0.001, alpha=0.2)
ggsave(p1, filename="images/UMAP.png", height=5, width=6)
# Test UMAP
DimPlot(dat_sub, reduction = "umap", group.by = "orig.clusters")

# Save 
dat_sc <- dat_sub
save(dat_sc, file="data/dat_singlecell.RData")
