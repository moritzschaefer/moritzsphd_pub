from gffutils.biopython_integration import to_seqfeature

from moritzsphd.data import grc_annotation, grc_seqs


def transcript_length(transcript, features=['exon']):
    db = grc_annotation('mesc')
    children = db.children(
        transcript, featuretype=features)
    return sum(child.end - child.start + 1 for child in children)


def transcript_sequence(transcript, features=['exon']):
    db = grc_annotation('mesc')
    seqs = grc_seqs('mesc')
    children = db.children(
        transcript, featuretype=features)
    return ''.join([str(to_seqfeature(child).extract(seqs[child.chrom]).seq.upper()) for child in children])
