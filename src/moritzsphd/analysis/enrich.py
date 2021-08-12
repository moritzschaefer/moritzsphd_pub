'''
simple functions to perform GO enrichment analysis

* TODO integrate GOMCL into my go analysis in moritzsphd
GOMCL: a toolkit to cluster, evaluate, and
extract non-redundant associations of Gene
Ontology-based functions
'''
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import gseapy as gp


def enrichr(gene_list, gene_sets=None, database=None, filter_padj=0.1, adjust_texts=False):
    '''
    Query the enrichR API with help of gseapy

    gene_list: A list of gene symbols(!)
    gene_sets: Provide a list of GO terms
    database: Alternatively use ALL terms from a database
              e.g. 'Mouse'
    '''

    if gene_sets is None and not database:
        gene_sets = ['WikiPathways_2019_Mouse', 'GO_Biological_Process_2018',
                     'GO_Molecular_Function_2018', 'GO_Cellular_Component_2018',
                     'KEGG_2019_Mouse', 'BioPlanet_2019', 'Reactome_2016',
                     'BioPlex_2017', 'CCLE_Proteomics_2020', 'Mouse_Gene_Atlas']
    elif database:
        gene_sets = gp.get_library_name(database=database)
    else:
        pass  # just use gene_sets

    enr = gp.enrichr(gene_list=gene_list,
                     gene_sets=gene_sets,
                     organism='mouse',
                     outdir='/tmp/',
                     cutoff=0.1  # test dataset, use lower value from range(0,1)
    )
    f, ax = plt.subplots(figsize=(13, 9))
    df = enr.results.loc[enr.results['Adjusted P-value'] < filter_padj].copy()
    df['-log10(adjusted p-value)'] = -np.log10(df['Adjusted P-value'])
    df['log2(#GO-target-genes)'] = df.Genes.apply(lambda genes: np.log2(len(genes.split(';'))))
    x = '-log10(adjusted p-value)'
    y = 'log2(#GO-target-genes)'
    sns.scatterplot(data=df, x=x, y=y, hue='Gene_set', s=150, alpha=0.8, ax=ax)
    texts = [ax.text(term[x], term[y], term['Term']) for i, term in df.iterrows()]
    if adjust_texts:
        from adjustText import adjust_text
        print(adjust_text(texts, arrowprops=dict(arrowstyle="->", color='r', lw=0.5)))

    return df, f
