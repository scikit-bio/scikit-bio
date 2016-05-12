# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.metadata._repr import _MetadataReprBuilder


class PairwiseDistancesReprBuilder(_MetadataReprBuilder):
    def _process_header(self):
        cls_name = self._obj.__class__.__name__
        self._lines.add_line(cls_name)
        self._lines.add_separator()

    def _process_data(self):
        ids = sorted(self._obj.ids)
        if len(ids) <= 5:
            self._lines.add_lines(self._format_id_counts(ids))
        else:
            self._lines.add_lines(self._format_id_counts(ids[:2]))
            self._lines.add_line('...')
            self._lines.add_lines(self._format_id_counts(ids[-2:]))

    def _format_id_counts(self, ids):
        lines = []
        for id1 in ids:
            count = 0
            for id2 in self._obj.ids:
                if id1 == id2:
                    continue
                if (id1, id2) in self._obj:
                    count += 1
            lines.append('%r: %d "between" distance%s' %
                         (id1, count, '' if count == 1 else 's'))
        return lines
