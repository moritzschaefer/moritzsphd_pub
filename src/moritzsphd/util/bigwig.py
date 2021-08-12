import re

import pandas as pd
import pyBigWig


# TODO I should calculate this on a per-transcript basis and then take the max- of all transcripts per gene
def tss_means(fn, genes, db, region='tss-3000+3000'):
    '''
    fn: bigwig-file
    genes: list of genes (from gffutils.FeatureDB)
    db: gffutils.FeatureDB
    region: either 'gene' or 'tss+<x>-<y>
    '''
    bw = pyBigWig.open(fn)

    def _gene_tss_mean(gene):
        if gene.chrom not in [str(v+1) for v in range(19)] + ['X', 'Y']:
            return None

        if region == 'gene':
            start, end = gene.start, gene.end  # tss is inaccurate, but I don't care for now..
        else:
            before, after = re.match('tss-([0-9]+)\+([0-9]+)', region).groups()
            if gene.strand == '+':
                start, end = gene.start-int(before), gene.start+int(after)
            else:
                start, end = gene.end-int(after), gene.end+int(before)
        try:
            # this computes the mean. As we are talking about log2ratios (bamCompare) here, this is not the best thing
            # to do, however we are stuck with it now
            return bw.stats(f'chr{gene.chrom}', max(1, start), min(end, bw.chroms(f'chr{gene.chrom}')))[0]
        except Exception:
            print(start, end, gene.chrom, fn)
            raise

    return pd.Series(data=[_gene_tss_mean(gene) for gene in genes], index=[gene.id for gene in genes])
