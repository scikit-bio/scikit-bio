# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.diversity._util import (
    _validate_counts_vector,
    _validate_taxa_and_tree,
    _vectorize_counts_and_tree,
    _check_taxa_alias,
)


def _setup_pd(counts, taxa, tree, validate, rooted, single_sample):
    if validate:
        if single_sample:
            # only validate count if operating in single sample mode, they
            # will have already been validated otherwise
            counts = _validate_counts_vector(counts)
            _validate_taxa_and_tree(counts, taxa, tree, rooted)
        else:
            _validate_taxa_and_tree(counts[0], taxa, tree, rooted)

    counts_by_node, _, branch_lengths = _vectorize_counts_and_tree(counts, taxa, tree)

    return counts_by_node, branch_lengths


def _faith_pd(counts_by_node, branch_lengths):
    """Calculate Faith's phylogenetic diversity (Faith's PD) metric.

    Parameters
    ----------
    counts_by_node : ndarray of shape (n_samples, n_nodes)
        Total counts/abundances of taxa descending from individual nodes of the tree.
    branch_lengths : ndarray of shape (n_nodes,)
        Branch lengths of corresponding nodes of the tree.

    Returns
    -------
    float
        Faith's phylogenetic diversity (PD).

    """
    return (branch_lengths * (counts_by_node > 0)).sum()


def faith_pd(counts, taxa=None, tree=None, validate=True, otu_ids=None):
    r"""Calculate Faith's phylogenetic diversity (Faith's PD) metric.

    The Faith's PD metric is defined as:

    .. math::

       PD = \sum_{b \in T \sqcup R} l(b)

    where :math:`T` is a minimum set of branches (:math:`b`) in a rooted tree
    that connect all taxa in a community. :math:`R` is a set of branches from
    the lowest common ancestor (LCA) of the taxa to the root of the tree.
    :math:`PD` is the sum of lengths (:math:`l`) of branches in both sets.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vectors of counts/abundances of taxa for one sample.
    taxa : list, np.array
        Vector of taxon IDs corresponding to tip names in ``tree``. Must be the
        same length as ``counts``. Required.
    tree : skbio.TreeNode
        Tree relating taxa. The set of tip names in the tree can be a superset
        of ``taxa``, but not a subset. Required.
    validate : bool, optional
        If ``False``, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results or error
        messages that are hard to interpret, so this step should not be
        bypassed if you're not certain that your input data are valid. See
        :mod:`skbio.diversity` for the description of what validation entails
        so you can determine if you can safely disable validation.
    otu_ids : list, np.array
        Alias of ``taxa`` for backward compatibility. Deprecated and to be
        removed in a future release.

    Returns
    -------
    float
        Faith's phylogenetic diversity (PD).

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails. Exact error will depend on what was invalid.

    See Also
    --------
    phydiv
    skbio.diversity
    skbio.diversity.alpha_diversity

    Notes
    -----
    Faith's phylogenetic diversity, often referred to as PD, was originally
    described in [1]_. It is the total phylogenetic branch length spanning
    all taxa of a community.

    It was clarified that the calculation should extend to the root of the
    tree [2]_, such that a single-taxon community will not have PD = 0. The
    root should be ancestral to all taxa being considered in the study, but
    does not need to be the origin of life. One should choose the root
    according to the scope of the study.

    Unrooted and abundance-weighted variants of PD are implemented in
    ``phydiv``.

    Several other metrics, such as evolutionary history (EH) [3]_ and
    functional diversity (FD) [4]_, are equivalent to PD in calculation.

    If computing Faith's PD for multiple samples, using
    ``skbio.diversity.alpha_diversity`` will be much faster than calling this
    function individually on each sample.

    This implementation of Faith's PD is based on the array-based
    implementation of UniFrac described in [5]_.

    This implementation differs from that in PyCogent (and therefore QIIME
    versions less than 2.0.0) by imposing a few additional restrictions on the
    inputs. First, the input tree must be rooted. In PyCogent, if an unrooted
    tree was provided that had a single trifurcating node (a newick convention
    for unrooted trees) that node was considered the root of the tree. Next,
    all taxa must be tips in the tree. PyCogent would silently ignore taxa that
    were not present the tree. To reproduce Faith PD results from PyCogent with
    scikit-bio, ensure that your PyCogent Faith PD calculations are performed
    on a rooted tree and that all taxa are present in the tree.

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    .. [2] Faith, D. P., & Baker, A. M. (2006). Phylogenetic diversity (PD)
       and biodiversity conservation: some bioinformatics challenges.
       Evolutionary bioinformatics, 2, 117693430600200007.

    .. [3] Nee, S., & May, R. M. (1997). Extinction and the loss of
       evolutionary history. Science, 278(5338), 692-694.

    .. [4] Petchey OL, Gaston KJ. Functional diversity (FD), species richness
       and community composition. Ecology letters. 2002 May;5(3):402-11.

    .. [5] Hamady M, Lozupone C, Knight R. Fast UniFrac: facilitating high-
       throughput phylogenetic analyses of microbial communities including
       analysis of pyrosequencing and PhyloChip data. ISME J. 4(1):17-27
       (2010).

    Examples
    --------
    Assume we have the following abundance data for a sample ``u``, represented
    as a counts vector. These counts represent the number of times specific
    taxa were observed in the sample.

    >>> u_counts = [1, 0, 0, 4, 1, 2, 3, 0]

    Because Faith PD is a phylogenetic diversity metric, we need to know which
    taxon each count corresponds to, which we'll provide as ``taxa``.

    >>> taxa = ['U1', 'U2', 'U3', 'U4', 'U5', 'U6', 'U7', 'U8']

    We also need a phylogenetic tree that relates the taxa to one another.

    >>> from io import StringIO
    >>> from skbio import TreeNode
    >>> tree = TreeNode.read(StringIO(
    ...                      '(((((U1:0.5,U2:0.5):0.5,U3:1.0):1.0):0.0,'
    ...                      '(U4:0.75,(U5:0.5,((U6:0.33,U7:0.62):0.5'
    ...                      ',U8:0.5):0.5):0.5):1.25):0.0)root;'))

    We can then compute the Faith PD of the sample.

    >>> from skbio.diversity.alpha import faith_pd
    >>> pd = faith_pd(u_counts, taxa, tree)
    >>> print(round(pd, 2))
    6.95

    """
    taxa = _check_taxa_alias(taxa, tree, otu_ids)

    counts_by_node, branch_lengths = _setup_pd(
        counts, taxa, tree, validate, rooted=True, single_sample=True
    )

    return _faith_pd(counts_by_node, branch_lengths)


def _phydiv(counts_by_node, branch_lengths, rooted, weight):
    """Calculate generalized phylogenetic diversity (PD) metrics.

    Parameters
    ----------
    counts_by_node : ndarray of shape (n_samples, n_nodes)
        Total counts/abundances of taxa descending from individual nodes of the tree.
    branch_lengths : ndarray of shape (n_nodes,)
        Branch lengths of corresponding nodes of the tree.
    rooted : bool
        Whether the metric is calculated considering the root of the tree.
    weight : bool or float
        Whether and to what degree branch lengths should be weighted by the relative
        abundance of taxa descending from the branch.

    Returns
    -------
    float
        Phylogenetic diversity (PD).

    """
    # select branches connecting taxa
    included = counts_by_node > 0

    # get total counts
    counts_sum = counts_by_node.max()
    if counts_sum == 0.0:
        return 0.0

    # in unrooted mode, remove branches to root
    if rooted is False:
        included &= counts_by_node < counts_sum

    # in unweighted mode, simply sum branch lengths
    if not weight:
        return (branch_lengths * included).sum()

    # get relative abundances
    fracs_by_node = counts_by_node / counts_sum

    # calculate balances in unrooted mode
    if rooted is False:
        fracs_by_node = 2 * np.minimum(fracs_by_node, 1 - fracs_by_node)

    # raise relative abundances to the power of theta
    if isinstance(weight, float) and weight < 1.0:
        fracs_by_node **= weight

    return (branch_lengths * fracs_by_node).sum()


def phydiv(
    counts, taxa=None, tree=None, rooted=None, weight=False, validate=True, otu_ids=None
):
    r"""Calculate generalized phylogenetic diversity (PD) metrics.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vectors of counts/abundances of taxa for one sample.
    taxa : list, np.array
        Vector of taxon IDs corresponding to tip names in ``tree``. Must be the
        same length as ``counts``. Required.
    tree : skbio.TreeNode
        Tree relating taxa. The set of tip names in the tree can be a superset
        of ``taxa``, but not a subset. Required.
    rooted : bool, optional
        Whether the metric is calculated considering the root of the tree. By
        default, this will be determined based on whether the input tree is
        rooted. However, one can override it by explicitly specifying ``True``
        (rooted) or ``False`` (unrooted).
    weight : bool or float, optional
        Whether branch lengths should be weighted by the relative abundance of
        taxa descending from the branch (default: ``False``). A float within
        [0, 1] indicates the degree of partial-weighting (0: unweighted, 1:
        fully-weighted).
    validate: bool, optional
        Whether validate the input data. See ``faith_pd`` for details.
    otu_ids : list, np.array
        Alias of ``taxa`` for backward compatibility. Deprecated and to be
        removed in a future release.

    Returns
    -------
    float
        Phylogenetic diversity (PD).

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails. Exact error will depend on what was invalid.

    See Also
    --------
    faith_pd
    skbio.diversity
    skbio.diversity.alpha_diversity

    Notes
    -----
    Phylogenetic diversity (PD) metrics measure the diversity of a community
    with consideration of the phylogenetic relationships among taxa. In
    general, PD is the sum of branch lengths spanning across taxa, optionally
    weighted by their abundance.

    The most widely-adopted PD metric, Faith's PD [1]_, is defined as:

    .. math::

       PD = \sum_{b \in T \sqcup R} l(b)

    where :math:`T` is a minimum set of branches (:math:`b`) in a rooted tree
    that connect all taxa in a community. :math:`R` is a set of branches from
    the lowest common ancestor (LCA) of the taxa to the root of the tree.
    :math:`PD` is the sum of lengths (:math:`l`) of branches in both sets.

    It is equivalent to ``pd(..., rooted=True, weight=False)``.

    A variant of PD, which does not include the root in the calculation, was
    referred to by some authors as unrooted phylogenetic diversity (uPD) [2]_,
    as in contrast to rooted phylogenetic diversity (rPD, i.e., Faith's PD).
    uPD is defined as:

    .. math::

       PD = \sum_{b \in T} l(b)

    It is equivalent to ``pd(..., rooted=False, weight=False)``.

    See ``faith_pd`` for a discussion of the root.

    PD (rooted or unrooted) considers only the presence of taxa. Therefore, it
    can be considered as the phylogenetic generalization of species richness.
    However, there are advantages of incorporating abundance information in the
    measurement [3]_. A generalized framework of abundance-weighted PD is
    provided in [4]_.

    Abundance-weighted rooted PD (equivalent to :math:`RBWPD_{1}` described in
    [4]_) is analogous to the :math:`PD_{aw}` metric originally described in
    [5]_ with a multiplier. It is defined as:

    .. math::

       PD = \sum_{b \in T \sqcup R} l(b) p(b)

    where :math:`p` is the sum of relative (proportional) abundances of taxa
    descending from branch (:math:`b`).

    It is equivalent to ``pd(..., rooted=True, weight=True)``.

    Abundance-weighted unrooted PD (equivalent to :math:`BWPD_{1}` described in
    [4]_) is analogous to the :math:`\delta nPD` metric originally described in
    [6]_ with a multiplier. It is defined as:

    .. math::

       PD = 2 \sum_{b \in T} l(b) \min(p(b),1-p(b))

    In which the term :math:`2\min(p(b),1-p(b))` is the lesser of the relative
    abundance of descending taxa on either side of a branch, multiplied by two.
    It is referred to as the "balance" of taxon abundance in [4]_.

    It is equivalent to ``pd(..., rooted=False, weight=True)``.

    The contribution of taxon abundance to the metric can be adjusted using the
    ``weighted`` parameter when it is a float within [0, 1]. This factor was
    referred to as :math:`\theta` in [4]_. The metric, :math:`BWPD_{\theta}`,
    referred to as the balance-weighted phylogenetic diversity in [4]_, is
    defined as:

    .. math::

       PD = \sum_{b \in T} l(b) (2\min(p(b),1-p(b)))^\theta

    It is equivalent to ``pd(..., rooted=False, weight=theta)``.

    This metric falls back to unweighted PD when :math:`\theta=0` or fully-
    weighted PD when :math:`\theta=1`. The original publication tested
    :math:`\theta=0.25` or :math:`0.5` [4]_.

    The parameter :math:`\theta` is analogous to the parameter :math:`\alpha`
    in the generalized UniFrac metric [7]_.

    Likewise, the rooted version of balance-weighted phylogenetic diversity,
    :math:`RBWPD_{\theta}` [4]_ (although "balance" is not involved), is
    defined as:

    .. math::

        PD = \sum_{b \in T \sqcup R} l(b) p(b)^\theta

    It is equivalent to ``pd(..., rooted=True, weight=theta)``.

    It is important to report which metric is used. For practical perspective,
    we recommend the following denotions:

    - :math:`rPD`: rooted, unweighted PD (Faith's PD [1]_).
    - :math:`uPD`: unrooted, unweighted PD (uPD [2]_).
    - :math:`rPD_{w}`: rooted, weighted PD (analogous to
      :math:`PD_{aw}` [3]_).
    - :math:`uPD_{w}`: unrooted, weighted PD (analogous to
      :math:`\delta nPD` [4]_).
    - :math:`rPD_{w\theta}`: rooted, weighted PD with parameter
      :math:`\theta` (:math:`RBWPD_{\theta}` [5]_).
    - :math:`uPD_{w\theta}`: unrooted, weighted PD with parameter
      :math:`\theta` (:math:`BWPD_{\theta}` [5]_).

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    .. [2] Pardi, F., & Goldman, N. (2007). Resource-aware taxon selection for
       maximizing phylogenetic diversity. Systematic biology, 56(3), 431-444.

    .. [3] Chao, A., Chiu, C. H., & Jost, L. (2016). Phylogenetic diversity
       measures and their decomposition: a framework based on Hill numbers.
       Biodiversity Conservation and Phylogenetic Systematics, 14, 141-72.

    .. [4] McCoy, C. O., & Matsen IV, F. A. (2013). Abundance-weighted
       phylogenetic diversity measures distinguish microbial community states
       and are robust to sampling depth. PeerJ, 1, e157.

    .. [5] Vellend, M., Cornwell, W. K., Magnuson-Ford, K., & Mooers, A. Ø.
       (2011). Measuring phylogenetic biodiversity. Biological diversity:
       frontiers in measurement and assessment, 194-207.

    .. [6] Barker, G. M. (2002). Phylogenetic diversity: a quantitative
       framework for measurement of priority and achievement in biodiversity
       conservation. Biological Journal of the Linnean Society, 76(2), 165-194.

    .. [7] Chen, J., Bittinger, K., Charlson, E. S., Hoffmann, C., Lewis, J.,
       Wu, G. D., ... & Li, H. (2012). Associating microbiome composition with
       environmental covariates using generalized UniFrac distances.
       Bioinformatics, 28(16), 2106-2113.

    """
    taxa = _check_taxa_alias(taxa, tree, otu_ids)

    # whether tree is rooted should not affect whether metric can be calculated
    # ; it is common unrooted PD is calculated on a rooted tree
    counts_by_node, branch_lengths = _setup_pd(
        counts, taxa, tree, validate, rooted=False, single_sample=True
    )

    # if not specified, determine whether metric should be calculated in rooted
    # mode according to the tree
    if rooted is None:
        rooted = len(tree.root().children) == 2

    # validate weight parameter
    if (
        not isinstance(weight, (bool, int, float))
        or (w_ := float(weight)) < 0.0
        or w_ > 1.0
    ):
        raise ValueError("Weight parameter must be boolean or within [0, 1].")

    return _phydiv(counts_by_node, branch_lengths, rooted, weight)
