'''Design oligos from sgRNA sequences'''

from skbio.sequence import DNA


def design_oligos(sgrna):
    '''
    >>> design_oligos('TCTTCAGACACTCCAGAAGA')
    ('CACCGTCTTCAGACACTCCAGAAGA', 'AAACTCTTCTGGAGTGTCTGAAGAC')
    '''
    fwd = 'CACCG' + sgrna
    rvs = 'AAAC' + str(DNA(sgrna).reverse_complement()) + 'C'

    return fwd, rvs
