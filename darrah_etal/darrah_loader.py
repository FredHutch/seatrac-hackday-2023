import numpy as np
import scanpy as sc
import pandas as pd

import sys
from os.path import join as opj

from fg_shared import _fg_data

data_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/darrah_etal')

def load_data():
    """Meta-data file contains cell-level data for Week 13 and Week 25"""
    md = pd.read_csv(opj(data_folder, 'updated_alexandria_metadata.txt'), sep='\t', low_memory=False)
    md = md.iloc[1:]

    """Clustering is provided separately for the Wk13 and Wk25 data"""
    clust = pd.read_csv(opj(data_folder, 'week13_clusters.txt'), sep='\t')
    clust = clust.iloc[1:]
    clust = clust.assign(X_13=clust['X'].astype(float),
                         Y_13=clust['Y'].astype(float))

    md = pd.merge(md, clust[['NAME', 'X_13', 'Y_13']], on='NAME', how='left')

    clust = pd.read_csv(opj(data_folder, 'week25_clusters.txt'), sep='\t')
    clust = clust.iloc[1:]
    clust = clust.assign(X_25=clust['X'].astype(float),
                         Y_25=clust['Y'].astype(float))

    """Merging here will create NA missing values for the cells that were not clustered amd thats expected
    (e.g., week25 clusters do not include week13 cells and the meta-data contains both cells)"""
    md = pd.merge(md, clust[['NAME', 'X_25', 'Y_25']], on='NAME', how='left')

    cts13 = pd.read_csv(opj(data_folder, 'Week13.Filtered.cells.txt.gz'), sep='\t')
    cts25 = pd.read_csv(opj(data_folder, 'Week25.Filtered.cells.txt.gz'), sep='\t')
    
    return md, cts13, cts25