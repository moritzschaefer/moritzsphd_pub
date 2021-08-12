'''
Functions to perform combined searches in various data sources
'''
from typing import List, Union

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from .data import (ensembl_release, gene_symbol_dict, geneid_to_name,
                   hek_mirna_expression, hek_mrna_expression,
                   kallisto_transcripts, mirtarbase, mrna_expression,
                   special_normalized_mirnas, tarbase, targetscan_all,
                   transcript_id_dict)
from .external.heap_ago2 import HeapAgo2
from .util import average_sample


def mrna_targeting_mirnas(
        gene_id,
        check_mirwalk=False,
        min_mirna_expression=30,
        max_ts_score=-0.35,
        filter_weak_and_clip=True,
        show_venn=False,
        cellline='mesc'
):
    '''
    Identify mRNAs targeting miRNAs

    :filter_weak_and_clip: if True, Weak and CLIP interactions are filtered out from MTB and TB
    '''
    mirnas = {}
    
    # alternative: targetscan_summary().query(f'`miRNA family` == "{mirna_seed}"')
    mirnas['targetscan'] = targetscan_all(cellline).loc[
        (targetscan_all(cellline).geneId == gene_id)
        & (targetscan_all(cellline)['context++ score'] < max_ts_score), 'miRNA']

    mtb_filter = (mirtarbase(cellline).geneId == gene_id)
    if filter_weak_and_clip:
        mtb_filter &= ~mirtarbase(cellline)['Support Type'].str.contains('Weak')

    mirnas['mirtarbase'] = mirtarbase(cellline).loc[mtb_filter, 'miRNA'].drop_duplicates()

    tb_filter = (tarbase(cellline).geneId == gene_id)
    if filter_weak_and_clip:
        tb_filter &= ~(tarbase(cellline).method.str.contains('CLIP'))
    mirnas['tarbase'] = tarbase(cellline).loc[tb_filter, 'mirna'].drop_duplicates()

    # AGO2 HEAP-CLIP
    heap_ago2_sites = HeapAgo2().mre_mirnas
    heap_mirnas = heap_ago2_sites.loc[heap_ago2_sites.gene_id == gene_id].mirna

    mirnas['heap_ago2'] = heap_mirnas

    if check_mirwalk:
        raise NotImplementedError
    if cellline == 'mesc':
        mirna_wt = average_sample(special_normalized_mirnas())['WT']
    else:
        mirna_wt = hek_mirna_expression()
    if min_mirna_expression > 0:
        for dbname in mirnas:
            mirnas[dbname] = mirnas[dbname][
                mirnas[dbname].map(mirna_wt.get) >= min_mirna_expression]

    if show_venn:
        from venn import venn
        plt.figure()
        venn({k: set(v) for k, v in mirnas.items()})
        plt.show()

    index = list(pd.concat(mirnas.values(), ignore_index=True))
    df = pd.DataFrame(index=index,
                      data={
                          k: [v in mirnas[k].tolist() for v in index]
                          for k in mirnas.keys()
                      })

    df['expression'] = mirna_wt.reindex(df.index)

    return df


def mirna_targets(mirna,
                  check_mirwalk=False,
                  min_target_expression=10,
                  max_ts_score=-0.35,
                  filter_weak_and_clip=True,
                  show_venn=False,
                  cellline='mesc'):
    '''
    mtb, tb, ts(, mirwalk?) targets expressed in mESCs
    :filter_weak_and_clip: if True, Weak and CLIP interactions are filtered out from MTB and TB
    '''
    targets = {}
    # alternative: targetscan_summary().query(f'`miRNA family` == "{mirna_seed}"')
    ts_targets = targetscan_all(cellline).loc[
        (targetscan_all(cellline).miRNA == mirna)
        & (targetscan_all(cellline)['context++ score'] < max_ts_score)]
    # convert gene symbol to geneid
    targets['targetscan'] = ts_targets['geneId']

    mtb_filter = (mirtarbase(cellline).miRNA == mirna)
    if filter_weak_and_clip:
        mtb_filter &= ~mirtarbase(cellline)['Support Type'].str.contains('Weak')
    targets['mirtarbase'] = mirtarbase(cellline).loc[mtb_filter, 'geneId'].drop_duplicates()

    tb_filter = (tarbase(cellline).mirna == mirna)
    if filter_weak_and_clip:
        tb_filter &= (tarbase(cellline).method.str.contains('CLIP'))
    targets['tarbase'] = tarbase(cellline).loc[tb_filter, 'geneId'].drop_duplicates()

    # AGO2 HEAP-CLIP
    heap_ago2_sites = HeapAgo2().mre_mirnas
    heap_mirnas = heap_ago2_sites.loc[heap_ago2_sites.mirna == mirna].gene_id
    
    targets['heap_ago2'] = heap_mirnas

    if check_mirwalk:
        raise NotImplementedError
    if cellline == 'mesc':
        mrna_wt = average_sample(mrna_expression()).set_index('Geneid')['WT']
    else:
        mrna_wt = hek_mrna_expression()
    if min_target_expression > 0:
        for dbname in targets:
            targets[dbname] = targets[dbname][
                targets[dbname].map(mrna_wt.get) >= min_target_expression]

    if show_venn:
        from venn import venn
        plt.figure()
        venn({k: set(v) for k, v in targets.items()})
        plt.show()
    index = list(set(pd.concat(targets.values(), ignore_index=True)))
    df = pd.DataFrame(index=index,
                      data={
                          k: [v in targets[k].tolist() for v in index]
                          for k in targets.keys()
                      })

    df['symbol'] = df.index.map(lambda g: geneid_to_name(g, cellline=cellline))

    df['expression'] = mrna_wt.reindex(df.index)

    return df


def gene_transcript_expression(genes: Union[str, List[str]], plot: bool = False):
    '''
    Return the expressed transcripts for the provided genes. Optionally plot an expression heatmap
    :genes: List of ensemble gene ids and/or gene symbols
    '''
    if type(genes) == str:
        genes = [genes]

    df = kallisto_transcripts().reset_index().copy()
    df['gene_id'] = df.apply(lambda row: transcript_id_dict().get(row['target_id'], None), axis=1)
    df.dropna(inplace=True)
    df['gene_name'] = df.apply(lambda row: ensembl_release.gene_name_of_gene_id(row.gene_id), axis=1)
    gene_ids = [gene_symbol_dict()[gene] if gene in gene_symbol_dict() else gene for gene in genes]

    df = df[df.gene_id.isin(gene_ids)]
    df.index = pd.MultiIndex.from_tuples(df.apply(lambda row: (row.gene_name, row.target_id), axis=1))
    df.drop(['target_id', 'gene_name', 'gene_id'], axis=1, inplace=True)
    if plot:
        plt.figure(figsize=(6, len(df) / 2.5))
        ax = sns.heatmap(df, annot=True, fmt='.2g', cbar_kws={'label': 'CPM'})
        ax.set_ylim(0, len(df))
    return df
