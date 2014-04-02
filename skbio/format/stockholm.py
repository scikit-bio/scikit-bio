#!/usr/bin/env python
r"""
Write stockholm formatted strings (:mod:`skbio.format.stockholm`)
=================================================================

.. currentmodule:: skbio.format.stockholm

This module provides functions for writing stockholm strings.

Constants
---------

.. autosummary::
   :toctree: generated/

    Stockholm


Functions
---------

.. autosummary::
   :toctree: generated/

    stockholm_from_alignment


"""


#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import namedtuple

Stockholm = namedtuple("Stockholm", ("aln", "GF", "GS", "GR", "GC"))
r"""
Stockholm named tuple for holding all information in a stockholm file

Attributes
----------
aln: Alignment object
    Holds the actual sequence alignment in the file
GF: dict
    Holds GF info in the format {feature: info}
GS: dict of dicts
    Holds GS info in the format {seqlabel: {feature: info}}
GR: dict of dicts
    Holds GR info in the format {seqlabel: {feature: info}}
GC: dict
    Holds GC info in the format {feature: info}

Notes
-----
Works with collectons.OrderedDict if the order of information matters
"""


def stockholm_from_alignment(stdata):
    r"""
    Parses a Stockholm named tuple into a string with stockholm format

    Parameters
    ----------
    stdata: Stockholm named tuple
        Stockholm formatted named tuple with alignment and all information

    Returns
    -------
    sto: string
        Stockholm formatted string containing all information passed

    Examples
    --------
    Basic Stockholm file:

    >>> from skbio.core.alignment import Alignment
    >>> from skbio.core.sequence import RNA
    >>> from skbio.format.stockholm import Stockholm, stockholm_from_alignment
    >>> seqs = [RNA("ACC--G-GGGU", identifier="seq1"),
    ...     RNA("TCC--G-GGGA", identifier="seq2")]
    >>> a1 = Alignment(seqs)
    >>> sto = Stockholm(a1, None, None, None, {"SS_cons": "(((....)))"})
    >>> stockholm_from_alignment(sto)
    '# STOCKHOLM 1.0\nseq1          ACC--G-GGGU\nseq2          TCC--G-GGGA\n
    #=GC SS_cons  (((....)))\n//'

    You can also use lists in features to have multi-line information:

    >>> from skbio.core.alignment import Alignment
    >>> from skbio.core.sequence import RNA
    >>> from skbio.format.stockholm import Stockholm, stockholm_from_alignment
    >>> seqs = [RNA("ACC--G-GGGU", identifier="seq1"),
    ...     RNA("TCC--G-GGGA", identifier="seq2")]
    >>> a1 = Alignment(seqs)
    >>> GF = {'BM': ['cmbuild  -F CM SEED', 'cmsearch  -Z 274931 -E 1000000']}
    >>> sto = Stockholm(a1, GR=None, GS=None, GF=GF,
    ...     GC={"SS_cons": "(((....)))"})
    >>> stockholm_from_alignment(sto)
    '# STOCKHOLM 1.0\n#=GF BM cmbuild  -F CM SEED\n
    #=GF BM cmsearch  -Z 274931 -E 1000000\nseq1          ACC--G-GGGU\n
    seq2          TCC--G-GGGA\n#=GC SS_cons  (((....)))\n//'

    However, lists of references and trees are treated differently:

    >>> from skbio.core.alignment import Alignment
    >>> from skbio.core.sequence import RNA
    >>> from skbio.format.stockholm import Stockholm, stockholm_from_alignment
    >>> seqs = [RNA("ACC--G-GGGU", identifier="seq1"),
    ...     RNA("TCC--G-GGGA", identifier="seq2")]
    >>> a1 = Alignment(seqs)
    >>> GF = {"RT": ["TITLE1",  "TITLE2"], "RL": ["J Mol Biol", "Cell"],
    ...       "RM": ["11469857", "12007400"], "NH": ["IMATREE", "IMATREETOO"],
    ...       "TN": ["Tree2", "Tree1"]}
    >>> sto = Stockholm(a1, GR=None, GS=None, GF=GF,
    ...     GC={"SS_cons": "(((....)))"})
    >>> stockholm_from_alignment(sto)
    '# STOCKHOLM 1.0\n#=GF RN [1]\n#=GF RM 11469857\n#=GF RT TITLE1\n
    #=GF RL J Mol Biol\n#=GF RN [2]\n#=GF RM 12007400\n#=GF RT TITLE2\n
    #=GF RL Cell\n#=GF TN Tree2\n#=GF NH IMATREE\n#=GF TN Tree1\n
    #=GF NH IMATREETOO\nseq1          ACC--G-GGGU\nseq2          TCC--G-GGGA\n
    #=GC SS_cons  (((....)))\n//'

    Notes
    -----
    If references are included in GF data, the RN lines are automatically
    generated if not provided.

    SPECIAL CASES:

    If there are multiple references, include information for each R* line as a
    list, with reference 0 information in position 0 for all lists, etc. This
    list will be broken up into the appropriate bits for each reference.

    If there are multiple trees included, use a list to store identifiers and
    trees, with position 0 holding identifier for tree in position 0, etc.
    """

    #find length of leader info needed to make file pretty
    #10 comes from the characters for '#=GF ' and the feature after label
    infolen = max(len(label) for label in stdata.aln.identifiers()) + 10

    GF_lines = []
    GS_lines = []
    GC_lines = []
    #NOTE: EVERYTHING MUST BE COERECED TO STR in case ints or floats are passed
    #add GF information if applicable
    if stdata.GF:
        for feature, value in stdata.GF.items():
            #list of features to skip and parse special later
            if feature in ("NH", "RC", "RM", "RN", "RA", "RL"):
                continue
            #list of features to parse special
            elif feature == "TN":
                #trees must be in proper order of identifier then tree
                ident = value if isinstance(value, list) else [value]
                tree = stdata.GF["NH"] if isinstance(stdata.GF["NH"], list) \
                    else [stdata.GF["NH"]]
                for ident, tree in zip(stdata.GF["TN"], stdata.GF["NH"]):
                    GF_lines.append(' '.join(["#=GF", "TN", str(ident)]))
                    GF_lines.append(' '.join(["#=GF", "NH", str(tree)]))
            elif feature == "RT":
                #make sure each reference block stays together
                #set up lists to zip in case some bits are missing
                #create rn list if needed: always have to have numbered refs
                rn = stdata.GF["RN"] if "RN" in stdata.GF else ["[%i]" % x for
                     x in range(1, len(value)+1)]
                rm = stdata.GF["RM"] if "RM" in stdata.GF else [0]*len(value)
                rt = stdata.GF["RT"] if "RT" in stdata.GF else [0]*len(value)
                ra = stdata.GF["RA"] if "RA" in stdata.GF else [0]*len(value)
                rl = stdata.GF["RL"] if "RL" in stdata.GF else [0]*len(value)
                rc = stdata.GF["RC"] if "RC" in stdata.GF else [0]*len(value)
                #order: RN, RM, RT, RA, RL, RC
                for n, m, t, a, l, c in zip(rn, rm, rt, ra, rl, rc):
                    GF_lines.append(' '.join(["#=GF", "RN", n]))
                    if m:
                        GF_lines.append(' '.join(["#=GF", "RM", str(m)]))
                    if t:
                        GF_lines.append(' '.join(["#=GF", "RT", str(t)]))
                    if a:
                        GF_lines.append(' '.join(["#=GF", "RA", str(a)]))
                    if l:
                        GF_lines.append(' '.join(["#=GF", "RL", str(l)]))
                    if c:
                        GF_lines.append(' '.join(["#=GF", "RC", str(c)]))
            else:
                #normal addition for everything else
                if not isinstance(value, list):
                    value = [value]
                for val in value:
                    GF_lines.append(' '.join(["#=GF", feature, str(val)]))

    #add GS information if applicable
    if stdata.GS:
        for seqname in stdata.GS:
            for feature in stdata.GS[seqname]:
                GS_lines.append(' '.join(["#=GS", seqname, feature,
                                         str(stdata.GS[seqname][feature])]))

    #add GC information if applicable
    if stdata.GC:
        for feature, value in stdata.GC.items():
            leaderinfo = ' '.join(["#=GC", feature])
            spacer = ' ' * (infolen - len(leaderinfo))
            GC_lines.append(spacer.join([leaderinfo, str(stdata.GC[feature])]))

    gr = stdata.GR if stdata.GR else {}

    sto_lines = ["# STOCKHOLM 1.0"] + GF_lines + GS_lines
    #create seq output along with GR info if applicable
    for label, seq in stdata.aln.iteritems():
        spacer = ' ' * (infolen - len(label))
        sto_lines.append(spacer.join([label, str(seq)]))
        #GR info added if exists for sequence
        if label in gr:
            for feature, value in gr[label].iteritems():
                leaderinfo = ' '.join(['#=GR', label, feature])
                spacer = ' ' * (infolen - len(leaderinfo))
                sto_lines.append(spacer.join([leaderinfo, value]))

    sto_lines.extend(GC_lines)
    #add final slashes to end of file
    sto_lines.append('//')

    return '\n'.join(sto_lines)
