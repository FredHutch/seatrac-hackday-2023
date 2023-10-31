import numpy as np
import scanpy as sc
import pandas as pd

import sys
from os.path import join as opj

from fg_shared import _fg_data

bigdata_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/bigdata')
wk4_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/gideon_etal/4week')
wk10_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/gideon_etal/10week')


def make_pseudo_bulk_4week():
    md = pd.read_csv(opj(wk4_folder, 'Updated4wk_alexandria_structured_metadata3.txt'), sep='\t')
    md = md.iloc[1:]

    bc = pd.read_csv(opj(wk4_folder, '4Week_barcodes.tsv'), sep='\t', header=None)
    feat = pd.read_csv(opj(wk4_folder, '4Week_features.tsv'), sep='\t', header=None)

    clust = pd.read_csv(opj(wk4_folder, '4Week_ClusteringDF.csv'))
    clust = clust.iloc[1:]
    clust = clust.assign(X=clust['X'].astype(float),
                         Y=clust['Y'].astype(float))

    md = pd.merge(md, clust, on='NAME', how='left')

    cts = sc.read_mtx(opj(bigdata_folder, 'gideon_4Week_countsmatrix.mtx'))

    out = {}
    for (did, ctype), gby in md.groupby(['donor_id', 'CellTypeAnnotations']):
         out[(did, ctype)] = np.asarray(cts[:, gby.index].X.sum(axis=1)).squeeze()
    out = pd.DataFrame(out, index=pd.Series(feat.values.squeeze(), name='gene'))
    out.columns = [f'D{i}_{j}' for i, j in out.columns]
    out.to_csv(opj(wk4_folder, 'pseudo_bulk_4week.csv'), index=True)


def make_pseudo_bulk_10week():
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

    """One cell per row, aligned with columns of cts1 and cts2 combined"""
    md = pd.merge(md, clust[['NAME', 'all_X', 'all_Y']], on='NAME', how='left')
    md = pd.merge(md, tclust[['NAME', 'T_X', 'T_Y']], on='NAME', how='left')

    
    cts1 = pd.read_csv(opj(bigdata_folder, 'gideon_10wk_counts_pt_1.csv.gz')).set_index('GENE')
    cts2 = pd.read_csv(opj(bigdata_folder, 'gideon_10wk_counts_pt_2.csv.gz')).set_index('GENE')
    
    """Genes (28155) x cells (109584)"""
    cts = pd.concat((cts1, cts2), axis=1)

    md = md.assign(sampleid=md.apply(lambda r: f'D{r["donor_id"]}_{r["biosample_id"]}', axis=1))

    out = {}
    for (ctype, sid), gby in md.groupby(['SpecificFinal', 'sampleid']):
         out[(ctype, sid)] = cts.iloc[:, gby.index].sum(axis=1)
    out = pd.DataFrame(out)
    out.columns = [f'{j}_{i}' for i,j in out.columns]
    out.to_csv(opj(wk10_folder, 'pseudo_bulk_10week.csv'), index=True)

if __name__ == '__main__':
    make_pseudo_bulk_4week()
    make_pseudo_bulk_10week()