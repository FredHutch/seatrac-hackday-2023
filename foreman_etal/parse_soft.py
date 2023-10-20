import pandas as pd

"""EXAMPLE CHUNK of SOFT file:
^SAMPLE = GSM7104701
!Sample_title = CD8_DFN0_19_AB2837
!Sample_geo_accession = GSM7104701
!Sample_status = Public on Mar 22 2023
!Sample_submission_date = Mar 19 2023
!Sample_last_update_date = Mar 22 2023
!Sample_type = SRA
!Sample_channel_count = 1
!Sample_source_name_ch1 = Lung
!Sample_organism_ch1 = Macaca mulatta
!Sample_taxid_ch1 = 9544
!Sample_characteristics_ch1 = tissue: Lung
!Sample_characteristics_ch1 = rna sample: AB2837
!Sample_characteristics_ch1 = cell type: CD8 granuloma T cell
!Sample_characteristics_ch1 = animal: DFN0
!Sample_characteristics_ch1 = granuloma: 19
!Sample_characteristics_ch1 = cfu: 2340
!Sample_characteristics_ch1 = cd8 t-cell_count: 15373
"""

def parse_file(fn):
    # fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/foreman_etal/GSE227653_family.soft')
    s = open(fn, 'r')

    keepers = ['Sample_title',
               'Sample_geo_accession',
               'Sample_source_name_ch1',
               'Sample_organism_ch1',
               'Sample_taxid_ch1',
               'Sample_characteristics_ch1']
    
    convert2int = ['granuloma', 'CFU', 'cd8 t-cell_count', 'cd4 t-cell_count']
    map_cols = {'rna sample':'biosample',
                'animal':'subjid',
                'cd8 t-cell_count':'cell_ct',
                'cd4 t-cell_count':'cell_ct'}

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
                    if key in convert2int:
                        value = int(value)
                else:
                    key = trimmed.replace('_ch1', '').replace('Sample_', '')

                key = map_cols.get(key, key)
                tmp.update({key : value})
    tab.append(tmp)
    out = pd.DataFrame(tab)
    return out


def load_data():
    data_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/foreman_etal/GSE227653_TPM_all.csv.gz')
    meta_fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/foreman_etal/GSE227653_family.soft.csv')

    meta = pd.read_csv(meta_fn)
    cts = pd.read_csv(data_fn)

    keep_cols = ['biosample', 'entrez_gene', 'gene_name_clc', 'Expression_value',
                  'TPM', 'RPKM', 'Exons', 'Gene_length', 'GeneID', 'Unique_gene_reads',
                  'Total_gene_reads', 'Transcripts_annotated', 'Detected_transcripts',
                  'Exon_length', 'Unique_exon_reads', 'Total_exon_reads']

    cts = pd.merge(cts[keep_cols], meta, how='left', on='biosample')
    return cts
    
if __name__ == '__main__':
    from fg_shared import _fg_data
    from os.path import join as opj

    fn = opj(_fg_data, 'SEATRAC/TB_hackday_2023/data/foreman_etal/GSE227653_family.soft')
    out = parse_file(fn)
    out.to_csv(fn + '.csv')
