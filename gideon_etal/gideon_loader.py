import numpy as np
import scanpy as sc
import pandas as pd

import sys
from os.path import join as opj

from fg_shared import _fg_data

wk4_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/gideon_etal/4week')
wk10_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/gideon_etal/10week')


def load_4week_data(load_raw=False):
    md = pd.read_csv(opj(wk4_folder, 'Updated4wk_alexandria_structured_metadata3.txt'), sep='\t')
    md = md.iloc[1:]

    bc = pd.read_csv(opj(wk4_folder, '4Week_barcodes.tsv'), sep='\t', header=None)
    feat = pd.read_csv(opj(wk4_folder, '4Week_features.tsv'), sep='\t', header=None)

    clust = pd.read_csv(opj(wk4_folder, '4Week_ClusteringDF.csv'))
    clust = clust.iloc[1:]
    clust = clust.assign(X=clust['X'].astype(float),
                         Y=clust['Y'].astype(float))

    md = pd.merge(md, clust, on='NAME', how='left')    

    if load_raw:
        cts = sc.read_mtx(opj(wk4_folder, '4Week_countsmatrix.mtx'))
    else:
        cts = sc.read_mtx(opj(wk4_folder, '4Week_datamatrix.mtx'))
    
    return md, bc, feat, cts

def load_10week_data(load_raw=False):
    md = pd.read_csv(opj(wk10_folder, 'Updated10wk_alexandria_structured_metadata10.txt'), sep='\t', low_memory=False)
    md = md.iloc[1:]

    clust = pd.read_csv(opj(wk10_folder, 'all_cells_umap.txt'), sep='\t')
    clust = clust.iloc[1:]
    clust = clust.assign(all_X=clust['X'].astype(float),
                         all_Y=clust['Y'].astype(float))

    tclust = pd.read_csv(opj(wk10_folder, 'T_Cells_UMAP.txt'), sep='\t')
    tclust = tclust.iloc[1:]
    tclust = tclust.assign(T_X=tclust['X'].astype(float),
                           T_Y=tclust['Y'].astype(float))

    md = pd.merge(md, clust[['NAME', 'all_X', 'all_Y']], on='NAME', how='left')
    md = pd.merge(md, tclust[['NAME', 'T_X', 'T_Y']], on='NAME', how='left')

    if load_raw:
        cts1 = pd.read_csv(opj(wk10_folder, 'counts_pt_1.csv.gz')).set_index('GENE')
        cts2 = pd.read_csv(opj(wk10_folder, 'counts_pt_2.csv.gz')).set_index('GENE')
    else:
        cts1 = pd.read_csv(opj(wk10_folder, 'lognormalized_pt_1.csv.gz')).set_index('GENE')
        cts2 = pd.read_csv(opj(wk10_folder, 'lognormalized_pt_2.csv.gz')).set_index('GENE')
    
    return md, pd.concat((cts1, cts2), axis=1)

