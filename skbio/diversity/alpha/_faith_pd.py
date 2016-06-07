# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.util._decorator import experimental
from skbio.diversity._util import (_validate_counts_vector,
                                   _validate_otu_ids_and_tree,
                                   _vectorize_counts_and_tree)


def _faith_pd(counts_by_node, branch_lengths):
    return (branch_lengths * (counts_by_node > 0)).sum()


@experimental(as_of="0.4.1")
def faith_pd(counts, otu_ids, tree, validate=True):
    """ Compute Faith's phylogenetic diversity metric (PD)

    Parameters
    ----------
    counts : 1-D array_like, int
        Vectors of counts/abundances of OTUs for one sample.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``counts``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results or error
        messages that are hard to interpret, so this step should not be
        bypassed if you're not certain that your input data are valid. See
        :mod:`skbio.diversity` for the description of what validation entails
        so you can determine if you can safely disable validation.

    Returns
    -------
    float
        The phylogenetic diversity (PD) of the samples.

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails. Exact error will depend on what was invalid.

    See Also
    --------
    skbio.diversity
    skbio.diversity.alpha_diversity

    Notes
    -----
    Faith's phylogenetic diversity, often referred to as PD, was originally
    described in [1]_.

    If computing Faith's PD for multiple samples, using
    ``skbio.diversity.alpha_diversity`` will be much faster than calling this
    function individually on each sample.

    This implementation differs from that in PyCogent (and therefore QIIME
    versions less than 2.0.0) by imposing a few additional restrictions on the
    inputs. First, the input tree must be rooted. In PyCogent, if an unrooted
    tree was provided that had a single trifurcating node (a newick convention
    for unrooted trees) that node was considered the root of the tree. Next,
    all OTU IDs must be tips in the tree. PyCogent would silently ignore OTU
    IDs that were not present the tree. To reproduce Faith PD results from
    PyCogent with scikit-bio, ensure that your PyCogent Faith PD calculations
    are performed on a rooted tree and that all OTU IDs are present in the
    tree.

    This implementation of Faith's PD is based on the array-based
    implementation of UniFrac described in [2]_.

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    .. [2] Hamady M, Lozupone C, Knight R. Fast UniFrac: facilitating high-
       throughput phylogenetic analyses of microbial communities including
       analysis of pyrosequencing and PhyloChip data.  ISME J. 4(1):17-27
       (2010).

    Examples
    --------
    Assume we have the following abundance data for a sample ``u``,
    represented as a counts vector. These counts represent the
    number of times specific Operational Taxonomic Units, or OTUs, were
    observed in the sample.

    >>> u_counts = [1, 0, 0, 4, 1, 2, 3, 0]

    Because Faith PD is a phylogenetic diversity metric, we need to know which
    OTU each count corresponds to, which we'll provide as ``otu_ids``.

    >>> otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5', 'OTU6', 'OTU7',
    ...            'OTU8']

    We also need a phylogenetic tree that relates the OTUs to one another.

    >>> from io import StringIO
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(StringIO(
    ...                      '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,'
    ...                      '(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62):0.5'
    ...                      ',OTU8:0.5):0.5):0.5):1.25):0.0)root;'))

    We can then compute the Faith PD of the sample.

    >>> from skbio.diversity.alpha import faith_pd
    >>> pd = faith_pd(u_counts, otu_ids, tree)
    >>> print(round(pd, 2))
    6.95

    """
    counts_by_node, branch_lengths = _setup_faith_pd(
        counts, otu_ids, tree, validate, single_sample=True)

    return _faith_pd(counts_by_node, branch_lengths)


def _setup_faith_pd(counts, otu_ids, tree, validate, single_sample):
    if validate:
        if single_sample:
            # only validate count if operating in single sample mode, they
            # will have already been validated otherwise
            counts = _validate_counts_vector(counts)
            _validate_otu_ids_and_tree(counts, otu_ids, tree)
        else:
            _validate_otu_ids_and_tree(counts[0], otu_ids, tree)

    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)

    return counts_by_node, branch_lengths
