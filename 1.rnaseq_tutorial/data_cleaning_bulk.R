library(tidyverse)
library(janitor)
library(edgeR)
library(limma)
library(RNAetc)
library(patchwork)

# Load Data
count <- read_csv("../darrah_etal/pseudo_bulk_wk13_wk25.csv") %>% 
  rename("gene"="...1")

# Select T-cells at week 25
count_tcell <- count %>% 
  select(gene, ends_with("T cell", ignore.case = FALSE)) %>% 
  select(gene, contains("WK25", ignore.case = FALSE))

# Create metadata from library names
meta <- data.frame(libID = colnames(count_tcell)[-1]) %>% 
  separate(libID, into=c("ptID","week","mtb","cell"), 
           sep="_", remove=FALSE) %>% 
  mutate(mtb = fct_recode(mtb, "Media"="STIMNO", "Mtb"="STIMYES"))

# Filter protein-coding genes
## Data base of genes
key <- read_tsv("data/Mmul_10_ensembl.txt") %>% 
  clean_names() %>% 
  filter(gene_type == "protein_coding") %>% 
  filter(gene_name %in% count_tcell$gene) %>% 
  select(gene_name, everything())

count_tcell_pc <- count_tcell %>% 
  filter(gene %in% key$gene_name)

# Create DGEList
## Order data
meta_ord <- meta %>% 
  arrange(libID)
count_tcell_pc_ord <- count_tcell_pc %>% 
  select(gene, all_of(meta_ord$libID)) %>% 
  arrange(gene) %>% 
  column_to_rownames("gene")
key_ord <- key %>% 
  filter(gene_name %in% count_tcell_pc$gene) %>% 
  arrange(gene_name) %>% 
  rename(gene=gene_name) %>% 
  #collapse multi annotations
  distinct() %>% 
  group_by(gene, gene_type) %>% 
  summarise(gene_stable_id = list(unique(gene_stable_id)))

dat <- DGEList(
  #count table
  counts=count_tcell_pc_ord,
  #metadata
  samples=meta_ord,
  #gene info
  genes=key_ord)

# Filter low abundance genes
## At least 0.1 CPM in at least 3 libraries
dat.abund <- RNAetc::filter_rare(dat, 
                                 min_CPM = 0.1, min.sample = 3,
                                 gene.var="gene")

## Visualize filter
p1 <- BIGpicture::plot_mv(dat, design = "~ mtb")
p2 <- BIGpicture::plot_mv(dat.abund, design = "~ mtb")
p1+p2

# Normalize
dat.abund.norm <- calcNormFactors(dat.abund, method = "TMM")
dat.abund.norm.voom <- voomWithQualityWeights(
  dat.abund.norm,
  design=model.matrix(~ mtb, 
                      data=dat.abund.norm$samples),
  plot=TRUE)

# Save
dat_tcell <- dat.abund.norm.voom
save(dat_tcell, file="data/dat_tcell.RData")
