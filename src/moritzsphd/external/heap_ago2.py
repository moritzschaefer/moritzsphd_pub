import re
from multiprocessing import Pool

import pandas as pd
# NOTE: grc_seqs uses GENCODE. I should probably switch to GRCm38_98 as on the server..
from moritzsphd.data import (ensembl_release, grc_annotation, grc_seqs,
                             mirbase_seqs, special_normalized_mirnas)
from moritzsphd.mirna_target_prediction import find_mres
from moritzsphd.util import average_sample
from moritzsphd.util.file import remote_file
from tqdm import tqdm
from tqdm.contrib.concurrent import process_map


def bed_export_mre_mirnas(filename, mre_mirnas):
    '''
    Compute MRE miRNAs and store them in a bed file
    '''
    try:
        import pybedtools
    except ImportError:
        print('Please install pybedtools')
        return
    # dont strip 'let-'
    if 'color' not in mre_mirnas.columns:
        mre_mirnas['color'] = '0,0,0'
    grouped_df_mres = mre_mirnas \
        .query('mre_type != "6mer"') \
        .groupby(['mre_type', 'start', 'chr', 'strand', 'offset', 'score', 'end']).agg({
            'mirna': lambda mres: 'miR-' + '/'.join([mre.lstrip('mmu-').lstrip('-miR') for mre in mres]),
            'color': 'first'}).reset_index()

    grouped_df_mres['label'] = grouped_df_mres.apply(lambda row: f'{row["mre_type"]}:{row["mirna"]}', axis=1)
    grouped_df_mres['score'] = grouped_df_mres['score'].astype('int').astype('str')
    grouped_df_mres['start'] = grouped_df_mres['start'].astype('int')
    grouped_df_mres['end'] = grouped_df_mres['end'].astype('int')
    grouped_df_mres['thickStart'] = grouped_df_mres['start']
    grouped_df_mres['thickEnd'] = grouped_df_mres['end']
    # grouped_df_mres['strand'] = '.'
    pybedtools.BedTool([v for v in grouped_df_mres[['chr', 'start', 'end', 'label', 'score', 'strand',
                                                    'thickStart', 'thickEnd', 'color']].values]) \
              .saveas(filename)

def mre_target(row):
    '''
    TODO: maybe check whether the strandness fits
    '''
    try:
        index, row = row
    except:
        pass
        # index = row.name

    # TODO simply add if hit['strand'] == row['strand']
    hits = grc_annotation().region(seqid=row['chr'].replace('chr', ''),
                                    start=row['start'], end=row['end'])
    hits = pd.DataFrame([{
        'gene_id': hit.attributes['gene_id'][0],
        'protein_coding': 'protein_coding' in hit.attributes['gene_biotype'],
        'feature_type': hit.featuretype
        } for hit in hits])
    if len(hits) == 0:
        return hits

    hits = hits[hits.protein_coding]  # ignore genes that are not protein_coding
    if len(hits) == 0:
        return hits
    hits = hits.groupby('gene_id').apply(lambda group: pd.Series({
        'is_cds': (group.feature_type == 'CDS').any(),
        'is_5putr': (group.feature_type == 'five_prime_utr').any(),
        'is_3putr': (group.feature_type == 'three_prime_utr').any(),
        # 'heap_mirna_index': index,
        **row
    }))
    if len(hits) == 0:
        return hits

    # drop intronic hits
    hits = hits.loc[hits.is_3putr | hits.is_5putr | hits.is_cds]
    return hits

def iterate_mirna_seqs(args):
    mirna_id, mirna_seq, row = args
    seq = grc_seqs()[str(row['seqnames']).replace('chr', '')][row['start']:row['end']].seq
    if row['strand'] == '-':
        seq = seq.reverse_complement()

    mre_6mer = str(mirna_seq[1:7].reverse_complement().back_transcribe())

    seed_matches = [m.start() for m in re.finditer(mre_6mer, str(seq))]

    # iterate over all identified 6mers
    mirnas = []
    for index in seed_matches:
        # see if the 6mer is a 8mer, 7merA1 or 7merm8
        try:
            a1 = seq[index + 6] == 'A'
        except IndexError:
            a1 = False
        try:
            m8 = str(seq[index -
                    1:index].reverse_complement()) == str(mirna_seq[7:8].back_transcribe())
        except IndexError:
            m8 = False

        if a1 and m8:
            mre_type = '8mer'
        elif a1:
            mre_type = '7merA1'
        elif m8:
            mre_type = '7merm8'
        else:  # discard 6mers later
            mre_type = '6mer'
        if row['strand'] == '-':
            end = row['start'] + (len(seq) - index) - 1 + int(m8)
            start = end - (6 + int(a1) + int(m8))
        else:
            start = row['start'] + index + 1 - int(m8)
            end = start + 6 + int(a1) + int(m8)
        mirnas.append({'mirna': mirna_id, 'mre_type': mre_type, 'chr': row['seqnames'].replace('chr', ''), 'start': start, 'end': end, 'strand': row['strand'], 'offset': index, 'score': row['score']})
    return mirnas


def peak_mirnas(row, positive=True, min_mirna_expression=0):
    if min_mirna_expression > 0:
        averaged_mirnas = average_sample(special_normalized_mirnas())['WT']
        if positive:
            mirnas = averaged_mirnas.index[averaged_mirnas > min_mirna_expression]
        else:
            mirnas = averaged_mirnas.index[averaged_mirnas < min_mirna_expression]
    else:
        mirnas = [mirna for mirna in mirbase_seqs().keys() if mirna.startswith('mmu-')]
    mirna_seqs = [(mirna, mirbase_seqs()[mirna], row) for mirna in mirnas]

    mres = pool.map(iterate_mirna_seqs, mirna_seqs)
    return [item for sublist in mres for item in sublist]


def missing_zic2_peak():
    '''
    manually add a peak for Zic2, because the transcript annotation is too
    short and misses the end of the 3'UTR
    [[file:~/wiki/gtd/reviews.org :ID:       4e989730-755c-4d9b-820b-f38651fe2f1f]]
    '''
    mirnas = peak_mirnas({'seqnames': '14', 'start': 122479890, 'end': 122479962, 'strand': '+', 'score': 5})
    df = pd.DataFrame(mirnas)
    df['gene_id'] = 'ENSMUSG00000061524'
    df['is_cds'] = False
    df['is_5putr'] = False
    df['is_3putr'] = True

    return df

pool = Pool(processes=1)  # TODO num processes should be defined wit num_threads as below..

class HeapAgo2:
    '''
    Computes mre_mirnas (access directly)
    '''
    def __init__(self, min_mirna_expression=10, mre_mirnas=None, num_threads=4):
        self.num_threads = num_threads
        filename = remote_file(
            'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139345&format=file&file=GSE139345_mESC_peaks.csv.gz',
            gunzip=True)
        self.df = pd.read_csv(filename, index_col=None)

        self.min_mirna_expression = min_mirna_expression

        tqdm.pandas()
        if mre_mirnas is None:
            self.mre_mirnas = self._gene_targets(self._mre_mirnas(positive=True))
            # glitch in zic2 annotation
            self.mre_mirnas = pd.concat((self.mre_mirnas, missing_zic2_peak())).reset_index(drop=True)
        else:
            self.mre_mirnas = mre_mirnas

    def _mre_mirnas(self, positive=True):
        '''
        Compute MRE miRNAs
        :positive: If false, use nonexpressed miRNAs!
        '''

        mre_mirnas = self.df.progress_apply(lambda row: peak_mirnas(row, positive), axis=1)
        df_mres = pd.DataFrame(list(mre_mirnas.explode().dropna()))

        return df_mres

    def _gene_targets(self, df):
        # df_mres = df.copy()
        # def _find_gene(row):
        #     try:
        #         return ensembl_release.genes_at_locus(contig=row['chr'][3:], position=row['start'], end=row['end'], strand=row['strand'])[0].gene_id
        #     except IndexError:
        #         return None

        # tmp = df.copy().reset_index()
        # tmp['chr'] = tmp['chr'].str.replace('chr', '')
        # bt = BedTool.from_dataframe(df[['chr', 'start', 'end', 'index']])

        # return pd.concat(df.progress_apply(_mirna_target, axis=1).tolist()) \
        results = process_map(mre_target, list(df.iterrows()), max_workers=self.num_threads, chunksize=100)
        results = [line for subdf in results for index, line in subdf.iterrows()]

        return pd.concat(results, axis=1, keys=[s.name for s in results], names=['gene_id']).T.reset_index()
