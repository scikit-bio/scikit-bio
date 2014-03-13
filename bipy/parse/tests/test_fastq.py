#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from cogent.util.unit_test import TestCase, main
from bipy.parse.fastq import MinimalFastqParser

data = {
"GAPC_0015:6:1:1259:10413#0/1":
        dict(seq='AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
             qual=r'````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'),
"GAPC_0015:6:1:1283:11957#0/1":
        dict(seq='TATGTATATATAACATATACATATATACATACATA',
             qual=r']KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb'),
"GAPC_0015:6:1:1284:10484#0/1":
        dict(seq='TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG',
             qual=r'UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_'),
"GAPC_0015:6:1:1287:17135#0/1":
        dict(seq='TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA',
             qual=r'^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K'),
"GAPC_0015:6:1:1293:3171#0/1":
        dict(seq="AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA",
             qual=r"b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c"),
"GAPC_0015:6:1:1297:10729#0/1":
        dict(seq="TAATGCCAAAGAAATATTTCCAAACTACATGCTTA",
             qual=r"T\ccLbb``bacc]_cacccccLccc\ccTccYL^"),
"GAPC_0015:6:1:1299:5940#0/1":
        dict(seq="AATCAAGAAATGAAGATTTATGTATGTGAAGAATA",
             qual=r"dcddbcfffdfffd`dd`^`c`Oc`Ybb`^eecde"),
"GAPC_0015:6:1:1308:6996#0/1":
        dict(seq="TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA",
             qual=r"a]`aLY`Y^^ccYa`^^TccK_X]\c\c`caTTTc"),
"GAPC_0015:6:1:1314:13295#0/1":
        dict(seq="AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT",
             qual=r"cLcc\\dddddaaYd`T```bLYT\`a```bZccc"),
"GAPC_0015:6:1:1317:3403#0/1":
        dict(seq="TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG",
             qual=r"\\\ZTYTSaLbb``\_UZ_bbcc`cc^[ac\a\T\ ".strip()) # had to add space
}
class ParseFastq(TestCase):
    def test_parse(self):
        """sequence and info objects should correctly match"""
        for label, seq, qual in MinimalFastqParser('data/fastq.txt'):
            self.assertTrue(label in data)
            self.assertEqual(seq, data[label]["seq"])
            self.assertEqual(qual, data[label]["qual"])
    

if __name__ == "__main__":
    main()
