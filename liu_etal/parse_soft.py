import pandas as pd

"""EXAMPLE CHUNK of SOFT file:
FROM GSE218157

^SAMPLE = GSM6735570
!Sample_title = 0NR IV pre
!Sample_geo_accession = GSM6735570
!Sample_status = Public on Jun 19 2023
!Sample_submission_date = Nov 16 2022
!Sample_last_update_date = Jun 19 2023
!Sample_type = SRA
!Sample_channel_count = 1
!Sample_source_name_ch1 = Whole Blood
!Sample_organism_ch1 = Macaca mulatta
!Sample_taxid_ch1 = 9544
!Sample_characteristics_ch1 = tissue: Whole Blood
!Sample_characteristics_ch1 = animal_id: 0NR
!Sample_characteristics_ch1 = pitt_id: 13417
!Sample_characteristics_ch1 = vax_group: IV
!Sample_characteristics_ch1 = vax_dose: 3.70E+07
!Sample_characteristics_ch1 = dose_group: high
!Sample_characteristics_ch1 = time after bcg: pre
!Sample_characteristics_ch1 = total_mtb_cfu: 30
!Sample_characteristics_ch1 = log_mtb_cfu: 1.491
!Sample_characteristics_ch1 = grans_nx: 1
!Sample_characteristics_ch1 = protect_outcome: protected
!Sample_molecule_ch1 = total RNA
!Sample_extract_protocol_ch1 = Whole blood was collected in PAXgene blood RNA tubes (Qiagen) and stored at -80C for batch processing at end of study. RNA was extracted using the PAXgene Blood RNA Tube kit (PreAnalytiX) as instructed. Globin mRNA was removed using GLOBINclear Kit (Life Technologies) and remaining mRNA concentration and quality was measured on an Agilent Bioanalyzer using an Agilent nano 6000 kit.
!Sample_extract_protocol_ch1 = Illumina-ready libraries were generated using NEBNext Ultra II RNA Preparation reagents (New England BioLabs).
!Sample_description = 0NR_65

FROM GSE218270

^SAMPLE = GSM6738257
!Sample_title = 13N022 high dose IV BCG pre
!Sample_geo_accession = GSM6738257
!Sample_status = Public on Jun 19 2023
!Sample_submission_date = Nov 17 2022
!Sample_last_update_date = Jun 19 2023
!Sample_type = SRA
!Sample_channel_count = 1
!Sample_source_name_ch1 = Whole Blood
!Sample_organism_ch1 = Macaca mulatta
!Sample_taxid_ch1 = 9544
!Sample_characteristics_ch1 = tissue: Whole Blood
!Sample_characteristics_ch1 = animal id: 13N022
!Sample_characteristics_ch1 = pitt id: 3218
!Sample_characteristics_ch1 = bcg dose_log10: 6.39
!Sample_characteristics_ch1 = dose group: high
!Sample_characteristics_ch1 = timeafterbcg: pre
!Sample_characteristics_ch1 = total mtb_cfu: 0
!Sample_characteristics_ch1 = log mtb_cfu: 0
!Sample_characteristics_ch1 = grans nx: 0
!Sample_characteristics_ch1 = log grans_nx: 0
!Sample_characteristics_ch1 = protect outcome: protected
!Sample_molecule_ch1 = total RNA
!Sample_extract_protocol_ch1 = Whole blood was collected in PAXgene blood RNA tubes (Qiagen) and stored at -80C for batch processing at end of study. RNA was extracted using the PAXgene Blood RNA Tube kit (PreAnalytiX) as instructed. Globin mRNA was removed using GLOBINclear Kit (Life Technologies) and remaining mRNA concentration and quality was measured on an Agilent Bioanalyzer using an Agilent nano 6000 kit.
!Sample_extract_protocol_ch1 = Illumina-ready libraries were generated using NEBNext Ultra II RNA Preparation reagents (New England BioLabs).
!Sample_description = 13N022_wkm4Pre

"""

def parse_file(fn):
    # fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218157_family.soft')
    # fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218270_family.soft')
    s = open(fn, 'r')

    keepers = ['Sample_title',
               'Sample_geo_accession',
               'Sample_source_name_ch1',
               'Sample_organism_ch1',
               'Sample_taxid_ch1',
               'Sample_characteristics_ch1',
               'Sample_description']

    convert2float = ['bcg dose_log10', 'total mtb_cfu', 'log mtb_cfu',
                     'total_mtb_cfu', 'log_mtb_cfu']

    map_cols = {'animal id':'subjid'}

    tab = []
    started = False
    for line in s.readlines():
        if line[:7] == '^SAMPLE':
            """Find each row that starts with "^SAMPLE" and then parse the next 17 rows if they match the
            keeper texts above"""
            if started:
                tab.append(tmp)
            started = True
            tmp = {}
            continue
        if started:
            trimmed = line[1:].split('=')[0].strip()
            if trimmed in keepers:
                value = line.split('=')[1].strip()
                if trimmed == 'Sample_characteristics_ch1':
                    key = value.split(':')[0].strip()
                    value = value.split(':')[1].strip()
                    if key in convert2float:
                        value = float(value)
                else:
                    key = trimmed.replace('_ch1', '').replace('Sample_', '')

                key = map_cols.get(key, key).replace(' ', '_')
                tmp.update({key : value})
    tab.append(tmp)
    out = pd.DataFrame(tab)
    return out

def load_data():
    """Load and merge meta-data into a pd DataFrame"""
    route_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218157_bcg_routecohort_processed.txt.gz')
    dose_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218270_ivbcg_dosecohort_processed.txt.gz')

    route_meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218157_family.soft.csv')
    dose_meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal/GSE218270_family.soft.csv')

    route = pd.read_csv(route_fn, sep='\t')
    dose = pd.read_csv(dose_fn, sep='\t')

    route_meta = pd.read_csv(route_meta_fn)
    dose_meta = pd.read_csv(dose_meta_fn)


    mmulatta_genes_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/seatrac-hackday-2023/mmulatta_to_human.tsv')
    mmulatta = pd.read_csv(mmulatta_genes_fn, sep='\t')


    route_ids = route.columns[1:]
    route_cts = pd.merge(route.rename({'Unnamed: 0': 'geneid'}, axis=1).set_index('geneid').stack().reset_index(),
                         mmulatta.rename({'Gene stable ID':'geneid'}, axis=1),
                         on='geneid', how='left')
    route_cts = route_cts.rename({'level_1':'sampleid',
                                  0:'ncounts'}, axis=1)
    route_cts = pd.merge(route_cts, route_meta.rename({'description':'sampleid'}, axis=1), on='sampleid', how='left')
    

    dose_ids = dose.columns[1:]
    dose_cts = pd.merge(dose.rename({'Unnamed: 0': 'geneid'}, axis=1).set_index('geneid').stack().reset_index(),
                         mmulatta.rename({'Gene stable ID':'geneid'}, axis=1),
                         on='geneid', how='left')
    dose_cts = dose_cts.rename({'level_1':'sampleid',
                                  0:'ncounts'}, axis=1)
    dose_cts = pd.merge(dose_cts, dose_meta.rename({'description':'sampleid'}, axis=1), on='sampleid', how='left')

    return route_cts, dose_cts
    
if __name__ == '__main__':
    from fg_shared import _fg_data
    from os.path import join as opj

    files = ['GSE218157_family.soft', 'GSE218270_family.soft']
    for fn in files:
        full_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/liu_etal', fn)
        out = parse_file(full_fn)
        out.to_csv(full_fn + '.csv')