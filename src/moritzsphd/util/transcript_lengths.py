import pandas as pd


# TODO here I could also use the transcript list from Daniel
def _filter_best_transcript(db, gene):
    transcripts = db.children(db[gene], featuretype='transcript')
    for appris_type in [
            'appris_principal_1',
            'appris_principal_2',
            'appris_principal_3',
            'appris_principal_4',
            'appris_principal_5',
            'appris_alternative_1',
            'appris_alternative_2',
    ]:
        for transcript in transcripts:
            try:
                if appris_type in transcript.attributes['tag']:
                    return transcript
            except KeyError:
                pass

    return next(db.children(db[gene], featuretype='transcript'))


def calculate_gene_lengths(genes, db=None):
    '''
    Calculate length of canonical mRNA
    '''
    if db == None:
        from moritzsphd.data import grc_annotation
        db = grc_annotation()
    gene_lengths = pd.Series(index=genes)
    for gene in genes:
        primary_transcript = _filter_best_transcript(db, gene)
        children = db.children(
            primary_transcript,
            featuretype=['exon'])
        gene_lengths.loc[gene] = sum(child.end - child.start + 1
                                     for child in children)

    return gene_lengths
