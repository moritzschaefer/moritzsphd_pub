'''
Dataset wrapper
# NOTE: grc_seqs and grc_annotation use GENCODE. I should probably switch to GRCm38_98 as on the server..
'''

import os
from functools import lru_cache, reduce

import pandas as pd
import pyensembl
from Bio import SeqIO

from moritzsphd.util import average_sample, generate_gffutils_db, remote_file

# TODO refactor everything with pyensembl!! use geneids WITHOUT version
ensembl_release = pyensembl.EnsemblRelease(98, pyensembl.species.mouse)


def cut_version(s: str):
    if type(s) is str:
        if '.' in s:
            return s[:s.find('.')]
        else:
            return s
    else:
        try:
            name = s.name
        except AttributeError:
            name = None
        return pd.Series([cut_version(v) for v in s], name=name)


@lru_cache()
def mrna_diffexp(group: str):
    '''
    diffexp of RNA-seq in our RNAi ggmutants

    group: defined in snakepipes
    '''

    group = group if '2i' in group else group.title()

    df = pd.read_csv(remote_file(
        f'cclab:/home/schamori/snakepipelines/rna-seq/DESeq2_{group}/DEseq_basic_DEresults_LFCshrunk.tsv'),
        # f'cclab:/home/schamori/moritzsphd/pipelines/rna-seq-star-deseq2/results/diffexp/{group}.diffexp.tsv'
                     sep='\t')
    return df


def geneid_to_name(gene_id, cellline='mesc'):
    '''
    Use like this:
    df['symbol'] = df['geneid'].apply(geneid_to_name)
    '''
    try:
        return ensembl_release.gene_name_of_gene_id(gene_id)
        # return grc_annotation(cellline)[gene_id].attributes['gene_name'][0]
    except (ValueError, KeyError) as e:
        if 'Run pyensembl install' in e.args[0]:
            raise
        else:
            return None
gid2n = geneid_to_name

def genename_to_id(gene_id, cellline='mesc'):
    '''
    Use like this:
    df['symbol'] = df['geneid'].apply(genename_to_id)
    '''
    try:
        return ensembl_release.gene_ids_of_gene_name(gene_id)[0]
        # return grc_annotation(cellline)[gene_id].attributes['gene_name'][0]
    except (ValueError, KeyError) as e:
        if 'Run pyensembl install' in e.args[0]:
            raise
        else:
            return None


gn2id = genename_to_id

@lru_cache()
def gene_symbol_dict(cellline='mesc'):
    print('gene_symbol_dict is deprecated. Use ensembl_release.gene_id_from_gene_name')
    d = {}
    for gene in ensembl_release.genes():
        d[gene.gene_name] = gene.gene_id

    return d


@lru_cache()
def transcript_id_dict():
    print('transcript_id_dict is deprecated. Use ensembl_release.gene_id_from_transcript_id')

    d = {}
    for transcript in ensembl_release.transcripts():
        d[transcript.transcript_id] = transcript.gene_id
        # assert f.id[:6] == 'ENSMUS',

    return d


@lru_cache()
def mrna_expression(): # normalization='.tpm'):
    '''
    '''
    df = pd.read_csv('/home/moritz/Projects/moritzsphd/pipelines/mesc-regulation/output/all.tpm.tsv',
        # remote_file(f'cclab:/home/schamori/moritzsphd/pipelines/rna-seq-star-deseq2/counts/all{normalization}.tsv'),
                     index_col=0,
                     sep='\t')
    # df.index = cut_version(df.index)
    # df = df.iloc[:, 5:]

    return df



@lru_cache()
def mirna_families(cellline='mesc'):
    download_prefix = 'mmu' if cellline == 'mesc' else 'vert'
    df = pd.read_csv(remote_file(
        f'http://www.targetscan.org/{download_prefix}_72/{download_prefix}_72_data_download/miR_Family_Info.txt.zip',
        'miR_Family_Info.txt'),
                     sep='\t')
    return df[df['Species ID'] == (10090 if cellline == 'mesc' else 9606)]  # mouse or human


@lru_cache()
def targetscan_summary():
    '''
    Per-family interactions
    '''
    df = pd.read_csv(remote_file(
        'http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.all_predictions.txt.zip',
        'Summary_Counts.txt'),
                     sep='\t')
    # df.drop(df.index[df['Cumulative weighted context++ score'] > -0.15],
    #         inplace=True)

    return df


@lru_cache()
def targetscan_conserved(cellline='mesc'):
    download_prefix = 'mmu' if cellline == 'mesc' else 'vert'
    df = pd.read_csv(remote_file(
        f'http://www.targetscan.org/{download_prefix}_72/{download_prefix}_72_data_download/Conserved_Site_Context_Scores.txt.zip',
        'Conserved_Site_Context_Scores.txt'),
                     sep='\t')
    # df.drop(df.index[df['Cumulative weighted context++ score'] > -0.15],
    #         inplace=True)

    df['geneId'] = cut_version(df['Gene ID'])
    df.drop(df.index[df['geneId'].isna()], inplace=True)

    return df.loc[df['Gene Tax ID'] == (10090 if cellline == 'mesc' else 9606), [
        'Gene Symbol', 'miRNA', 'context++ score',
        'context++ score percentile', 'weighted context++ score',
        'weighted context++ score percentile', 'geneId'
    ]]


@lru_cache()
def targetscan_nonconserved(cellline='mesc'):
    download_prefix = 'mmu' if cellline == 'mesc' else 'vert'
    df = pd.read_csv(remote_file(
        f'http://www.targetscan.org/{download_prefix}_72/{download_prefix}_72_data_download/Nonconserved_Site_Context_Scores.txt.zip',
        'Nonconserved_Site_Context_Scores.txt'),
                     sep='\t')
    # df.drop(df.index[df['Cumulative weighted context++ score'] > -0.15],
    #         inplace=True)

    df['geneId'] = cut_version(df['Gene ID'])
    # df['geneId'] = df['Gene Symbol'].apply(lambda gs: gene_symbol_dict(cellline).get(gs, None))
    df.drop(df.index[df['geneId'].isna()], inplace=True)

    df = df.loc[df['Gene Tax ID'] == (10090 if cellline == 'mesc' else 9606), [
        'Gene Symbol', 'miRNA', 'context++ score',
        'context++ score percentile', 'weighted context++ score',
        'weighted context++ score percentile', 'geneId'
    ]]
    # approximately 10% don't have a context++ score. delete them
    return df.loc[~df['context++ score'].isna()]


# TODO check if mutation of return value affects cache..
def targetscan_all(cellline='mesc'):
    df_cons = targetscan_conserved(cellline)
    df_noncons = targetscan_nonconserved(cellline)

    df = pd.concat([df_cons, df_noncons], ignore_index=True)
    df['conserved'] = False
    df.iloc[:len(df_cons), -1] = True

    return df


@lru_cache()
def special_normalized_mirnas():
    '''
    MiRNAs in our cell lines, normalized by tRNA, snRNA and snoRNA
    '''
    df = pd.read_csv(remote_file(
        'cclab:/home/schamori/moritzsphd/pipelines/srna-seq/deseq2/mirbase_mmu21.special_normalized_counts.tsv'
    ),
                     sep=' ')

    return df

@lru_cache()
def mirna_diffexp(group: str):
    '''
    diffexp of RNA-seq in our RNAi ggmutants

    group: defined in snakepipes
    '''

    df = pd.read_csv(remote_file(
        f'cclab:/home/schamori/moritzsphd/pipelines/srna-seq/results/diffexp.mirbase_mmu21/{group.lower()}.diffexp.tsv'
    ), sep=' ')

    return df


@lru_cache()
def grc_seqs(cellline='mesc'):
    '''
    Dict of strings for the chromosomes. Keys are 'chr1' to 'chrY'
    '''
    if cellline == 'mesc':
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz'
    else:
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh38.primary_assembly.genome.fa.gz'

    grc_file = remote_file(url, gunzip=True)
    fasta_sequences = SeqIO.parse(grc_file, 'fasta')
    return dict((key.replace('chr', ''), value) for key, value in SeqIO.to_dict(fasta_sequences).items())


@lru_cache()
def grc_annotation(cellline='mesc'):
    '''
    '''
    if cellline == 'mesc':
        # url = 'cclab:~/data/snakepipes/GRCm38_98/annotation/genes.gtf'
        url = 'ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz'
    else:
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gff3.gz'

    gff_filename = remote_file(url, gunzip=True)

    return generate_gffutils_db(gff_filename)


@lru_cache()
def mirbase_seqs(organism=''):
    mirna_seqs = SeqIO.parse(str(
        remote_file('ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz',
                    gunzip=True)), 'fasta')

    return {seq.id: seq.seq for seq in mirna_seqs if seq.id.startswith(organism)}



@lru_cache()
def primary_transcripts():
    db = grc_annotation('mesc')
    primary_transcripts = {}
    for gene in db.features_of_type('gene'):
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
                        primary_transcripts[gene.id] = transcript.id
                except KeyError:
                    pass
            else:
                continue
            break

    s = pd.Series(primary_transcripts)
    s.index.name = 'Geneid'
    return s


