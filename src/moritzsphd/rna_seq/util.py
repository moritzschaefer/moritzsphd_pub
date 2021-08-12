import numpy as np

PRIOR_COUNT = 0.25


def log_transform(self, counts, prior_count=PRIOR_COUNT):
    """Compute the log of the counts"""
    counts[counts == 0] = prior_count
    return np.log(counts)


def tpm(
        gene_counts,
        gene_lengths,
        transform_to_log=False,
        mean_fragment_length=None,
):
    """Normalize the DGEList to transcripts per million.
    Adapted from Wagner, et al. 'Measurement of mRNA abundance using RNA-seq data:
    RPKM measure is inconsistent among samples.' doi:10.1007/s12064-012-0162-3
    Read counts :math:`X_i` (for each gene :math:`i` with gene length :math:`\widetilde{l_j}` )
    are normalized as follows:
    .. math::
        TPM_i = \\frac{X_i}{\\widetilde{l_i}}\cdot \\
        \\left(\\frac{1}{\sum_j \\frac{X_j}{\widetilde{l_j}}}\\right) \cdot 10^6
    Args:
        gene_lengths: 1D array of gene lengths for each gene.
        transform_to_log: store log outputs
        prior_count:
        mean_fragment_length: mean fragment length of reads
            (optional)
    """

    # compute effective length not allowing negative lengths
    if mean_fragment_length:
        effective_lengths = (gene_lengths - mean_fragment_length).clip(min=1)
    else:
        effective_lengths = gene_lengths

    # how many counts per base
    base_counts = gene_counts / effective_lengths

    counts = 1e6 * base_counts / base_counts.sum()
    if transform_to_log:
        counts = log_transform(counts)

    return counts
