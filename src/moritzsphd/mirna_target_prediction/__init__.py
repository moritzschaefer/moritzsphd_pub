import re

import pandas as pd

from moritzsphd.data import mirbase_seqs


def find_mres(mrna_seq, mirnas):
    '''
    Returns a dataframe with binding sites for the given combination
    :mrna_seq: str of regulatory sequence
    :mirnas: mirna identifier (e.g. hsa-miR-290a) or list of mirna identifiers
    '''

    mrna_seq = mrna_seq.upper().replace('T', 'U')

    if not isinstance(mirnas, list):
        mirnas = [mirnas]

    if isinstance(mirnas[0], str):
        mirnas = [(mirna, mirbase_seqs()[mirna]) for mirna in mirnas]
    else:
        mirnas = [(seq.id, seq.seq) for seq in mirnas]

    rows = []
    for mirna, mirna_seq in mirnas:

        mirna_seeds = {
            mirna_seq[1:8].reverse_complement(): 'm8mer',
            mirna_seq[1:7].reverse_complement() + 'A': 'a1mer',
            mirna_seq[1:8].reverse_complement() + 'A': '8mer',
            mirna_seq[1:7].reverse_complement(): '6mer',
        }
        pattern = f'{mirna_seq[7:8].reverse_complement()}?{mirna_seq[1:7].reverse_complement()}A?'  # first and last characters are optional

        for result in re.finditer(pattern, mrna_seq):
            rows.append({
                'type': mirna_seeds[result.group(0)],
                'pos': result.start(),
                'mirna': mirna,
            })

    return pd.DataFrame(rows)
