"""
Reproducing clusters as found on mirbase. Use genome specific files like this
one: ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3
"""

import re
from urllib.request import urlretrieve

from ...config import DATA_DIR
from ...util import remote_file

MIRBASE_GFF = remote_file('ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3')
# '/home/moritz/data/mmu.gff3'


def read_mirna_gff(filename):
    """Read a gff3 file from mirbase and return a preprocessed dataframe

    Args:
        filename: Path to the gff3 file containing metadata for miRNAs.
    Returns:
        pandas.DataFrame: A dataframe containing all fields required by the
        clustering functions (name, seqname and start)
    """
    from pybedtools import BedTool
    mirna_df = BedTool(filename).to_dataframe().dropna()  # dropna removes header rows
    mirna_df['name'] = mirna_df.apply(
        lambda v: dict(s.split('=') for s in v.attributes.split(';'))['Name'],
        axis=1)
    # there are duplicates which we need to remove # TODO let-7k and let-7j still have duplicates...
    mirna_df.drop(
        mirna_df.index[mirna_df.apply(
            lambda v: dict(s.split('=') for s in v.attributes.split(';'))['ID'],
            axis=1).str.contains('_')],
        inplace=True)
    mirna_df['ID'] = mirna_df.apply(
        lambda v: dict(s.split('=') for s in v.attributes.split(';'))['ID'],
        axis=1)
    return mirna_df


def _find_neighbors(mirna_df, seqname, start, max_dist):
    neighbors = mirna_df[(mirna_df.seqname == seqname)
                         & (mirna_df.start > start - max_dist) &
                         (mirna_df.start < start + max_dist)]
    return set(neighbors['name'])


def _add_to_clusters(clusters, mirs):
    for i, c in enumerate(clusters):
        for mir in mirs:
            if mir in c:
                # add mirs to cluster
                clusters[i] = c | mirs
                return

    # if cluster was not found, create new cluster
    clusters.append(mirs)


def cluster_to_name(cluster, mirna_df=None):
    """Take a set of miRNAs and create a name from it

    Todos:
        - order by position
        - improve check for unique_names when manually naming clusters

    Args:
        cluster: set<string>. A set of strings (miRNA names)
        mirna_df: pd.Dataframe. If provided takes the spatially first and
            last mirnas for the name.

    Returns:
        string: A name miRNA cluster name
    """
    unique_names = set([
        re.match('mmu-(mir|miR|let)-(\d+\w?(-\d+)?)(-[35]p)?$', m).groups()[1]
        for m in cluster if '.' not in m
    ])

    if '1193' in unique_names:
        return 'miR-379~Mirg~3072'

    if '466a' in unique_names:
        return 'miR-466*-669*'

    if '337' in unique_names:
        return 'miR-337-136'

    if '293' in unique_names:
        return 'miR-290-295'

    if '183' in unique_names:
        return 'miR-183-182'

    if '17' in unique_names:
        return 'miR-17-92a'

    # TODO the sorting is bad..
    cluster_name = ','.join(sorted(unique_names)[:5])

    return f'miR-{cluster_name}'


def cluster_mirnas(mirna_data, max_dist=10000):
    """Cluster miRNAs based on proximity

    Args:
        mirna_data: Either a pandas DataFrame with columns (name, seqname
            and start) or a path(string) to a gff3 file containing such data
            (downloadable from mirBase)
        max_dist: Maximal distance for miRNAs to be considered in the same
            cluster. 10000 is the default and used by mirbase as well
    Returns
        list<set>: List of sets of mirnas that are close to each other
    """

    if isinstance(mirna_data, str):
        mirna_df = read_mirna_gff(mirna_data)
    else:
        mirna_df = mirna_data

    clusters = []
    for i, row in mirna_df.iterrows():
        # find a set(!) of neighbors (including oneself)
        neighbors = _find_neighbors(mirna_df, row.seqname, row.start, max_dist)
        _add_to_clusters(clusters, neighbors)

    return clusters


def mirbase_clusters():
    """Helper function which downloads mirbase and computes clusters"""
    mirna_df = read_mirna_gff(MIRBASE_GFF)
    clusters = cluster_mirnas(mirna_df)

    return clusters
