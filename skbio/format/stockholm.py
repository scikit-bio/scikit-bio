#!/usr/bin/env python
"""Writer for stockholm sequence format"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import namedtuple


"""
Stockholm named tuple takes the following format
aln   Alignment object
GS    {seqlabel: {feature: info}}
GF    {feature: info}
GR    {seqlabel: {feature: info}}
GC    {feature: info}
"""
Stockholm = namedtuple('Stockholm', ('aln', 'GS', 'GF', 'GR', 'GC'))


def stockholm_from_alignment(stdata):
    #find length of leader info needed to make file pretty
    #9 comes from the characters for '#=GF ' and the feature info after label
    infolen = max(len(label) for label in stdata.aln.identifiers()) + 10

    GF_lines = []
    GS_lines = []
    GC_lines = []

    #add GF information if applicable
    if stdata.GF:
        GF_lines = [' '.join(["#=GF", feature, str(value)])
                    for feature, value in stdata.GF.iteritems()]

    #add GS information if applicable
    if stdata.GS:
        for seqname in stdata.GS:
            for feature in stdata.GS[seqname]:
                GS_lines.append(' '.join(["#=GS", seqname, feature,
                                         str(stdata.GS[seqname][feature])]))

    #add GC information if applicable
    if stdata.GC:
        GC_lines = [' '.join(["#=GC", feature, str(value)])
                    for feature, value in stdata.GC.iteritems()]

    GR = stdata.GR if stdata.GR else {}

    sto_lines = ["# STOCKHOLM 1.0"] + GF_lines + GS_lines
    #create seq output along with GR info if applicable
    for label, seq in stdata.aln.iteritems():
        spacer = ' ' * (infolen - len(label))
        sto_lines.append(spacer.join([label, str(seq)]))
        #GR info created if exists for sequence
        if label in GR:
            for feature, value in GR[label].iteritems():
                leaderinfo = ' '.join(['#=GR', label, feature])
                spacer = ' ' * (infolen - len(leaderinfo))
                sto_lines.append(spacer.join([leaderinfo, value]))

    sto_lines.extend(GC_lines)
    #add final slashes to end of file
    sto_lines.append('//')

    return '\n'.join(sto_lines)
