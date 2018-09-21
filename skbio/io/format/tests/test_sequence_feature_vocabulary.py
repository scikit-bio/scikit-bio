# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.io.format._sequence_feature_vocabulary import (
    _parse_loc_str, _parse_section_default, _serialize_location)
from skbio.io import FileFormatError
from skbio.metadata import IntervalMetadata


class Tests(TestCase):
    def test_parse_section_default(self):
        lines = [
            ['FOO  blah blah',
             '     blah'],
            ['FOO=blah',
             '    blah'],
            ['FOO']]
        kwargs = [{'join_delimiter': '=', 'return_label': False},
                  {'label_delimiter': '=', 'join_delimiter': '',
                   'return_label': True},
                  {'label_delimiter': '=', 'join_delimiter': '=',
                   'return_label': True}]
        expects = ['blah blah=blah',
                   ('FOO', 'blahblah'),
                   ('FOO', '')]
        for i, j, k in zip(lines, kwargs, expects):
            self.assertEqual(k, _parse_section_default(i, **j))

    def test_parse_loc_str(self):
        examples = [
            '9',  # a single base in the presented sequence
            '3..8',
            '<3..8',
            '3..>8',
            'complement(3..>8)',
            'complement(join(3..>5,<7..9))',
            'join(J00194.1:1..9,3..8)',
            'join(3..8,J00194.1:1..9)',
            '1.9',
            '1^2']

        expects = [
            ([(8, 9)], [(False, False)], {'strand': '+'}),
            ([(2, 8)], [(False, False)], {'strand': '+'}),
            ([(2, 8)], [(True, False)],  {'strand': '+'}),
            ([(2, 8)], [(False, True)],  {'strand': '+'}),
            ([(2, 8)], [(False, True)],  {'strand': '-'}),
            ([(2, 5), (6, 9)], [(False, True), (True, False)],
             {'strand': '-'}),
            ([(2, 8)], [(False, False)], {'strand': '+'}),
            ([(2, 8)], [(False, False)], {'strand': '+'}),
            ([(0, 9)], [(False, False)], {'strand': '+'}),
            ([(0, 1)], [(False, False)], {'strand': '+'})]

        for example, expect in zip(examples, expects):
            parsed = _parse_loc_str(example)
            self.assertEqual(parsed, expect)

    def test_parse_loc_str_invalid(self):
        examples = [
            'abc',
            '3-8']
        for example in examples:
            with self.assertRaisesRegex(FileFormatError,
                                        r'Could not parse location string: '
                                        '"%s"' % example):
                _parse_loc_str(example)

    def test_serialize_location(self):
        imd = IntervalMetadata(9)
        i1 = imd.add([(0, 1)])
        self.assertEqual(_serialize_location(i1), '1')

        i2 = imd.add([(0, 2)], [(True, True)])
        self.assertEqual(_serialize_location(i2), '<1..>2')

        i3 = imd.add([(0, 2)], [(False, True)])
        self.assertEqual(_serialize_location(i3), '1..>2')

        i4 = imd.add([(0, 2)], [(True, False)])
        self.assertEqual(_serialize_location(i4), '<1..2')

        i5 = imd.add([(0, 2), (3, 9)], metadata={'strand': '-'})
        self.assertEqual(_serialize_location(i5),
                         'complement(join(1..2,4..9))')

        i6 = imd.add([(0, 2), (3, 9)],
                     [(True, False), (False, True)],
                     metadata={'strand': '-'})
        self.assertEqual(_serialize_location(i6),
                         'complement(join(<1..2,4..>9))')


if __name__ == '__main__':
    main()
