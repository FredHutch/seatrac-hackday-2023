import numpy as np
import pandas as pd
import sys
from os.path import join as opj

from fg_shared import _fg_data

data_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/darrah_etal')
bigdata_folder = opj(_fg_data, 'SEATRAC/TB_hackday_2023/bigdata')

def load_raw_data():
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

    cts13 = pd.read_csv(opj(bigdata_folder, 'darrah_Week13.Filtered.cells.txt'), sep='\t')
    cts25 = pd.read_csv(opj(bigdata_folder, 'darrah_Week25.Filtered.cells.txt'), sep='\t')
    
    """NOTE: the meta-data and cts data are not aligned by cellID so a MERGE is neccessary on md['NAME'] and cts columns"""
    return md, cts13, cts25

def load_pseudo_bulk():
    md = pd.read_csv(opj(data_folder, 'updated_alexandria_metadata.txt'), sep='\t', low_memory=False)
    md = md.iloc[1:]

    md_cols = ['sampleid', 'organ',
               'organ__ontology_label', 'vaccination', 'vaccination__ontology_label',
               'sex', 'is_living', 'sample_type', 'biosample_id', 'end_bias',
               'cell_type', 'cell_type__ontology_label', 'disease',
               'disease__ontology_label', 'donor_id', 'Stimulated', 'VaccineRoute',
               'VaccineRouteUnique', 'VaccineRouteUniqueGroup',
               'sequencing_instrument_manufacturer_model',
               'sequencing_instrument_manufacturer_model__ontology_label',
               'paired_ends', 'read_length', 'vaccination_route',
               'vaccination__time_since', 'vaccination__time_since__unit',
               'vaccination__time_since__unit_label', 'vaccination__dosage']
    md = md.assign(sampleid=md.apply(lambda r: f'D{r["donor_id"]}_WK{r["vaccination__time_since"]}_STIM{r["Stimulated"]}', axis=1))
    md[md_cols].drop_duplicates()

    cts = pd.read_csv(opj(data_folder, 'pseudo_bulk_wk13_wk25.csv'))
    return md, cts



