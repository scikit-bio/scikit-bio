#!/usr/bin/env python
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import csv
from itertools import izip
from StringIO import StringIO
import numpy as np

class MissingHeaderError(Exception):
    pass

class DistanceMatrix(object):

    @classmethod
    def from_file(cls, dm_f, delimiter='\t'):
        sample_ids = cls._extract_sample_ids(dm_f, delimiter)
        data = np.loadtxt(dm_f, delimiter=delimiter, skiprows=1,
                          usecols=range(1, len(sample_ids) + 1), ndmin=2)
        return cls(data, sample_ids)

    def __init__(self, data, sample_ids):
        data = np.asarray(data)
        sample_ids = tuple(sample_ids)
        self._validate(data, sample_ids)

        self.data = data
        self.sample_ids = sample_ids
        self._sample_index = self._index_list(self.sample_ids)

    @property
    def dtype(self):
        return self.data.dtype

    @property
    def shape(self):
        return self.data.shape

    def __str__(self):
        return '%dx%d distance matrix\nSample IDs: %s\n' % (self.shape[0],
                self.shape[1], ', '.join(self.sample_ids)) + str(self.data)

    def __eq__(self, other):
        eq = True

        # The order these checks are performed in is important to be as
        # efficient as possible. The check for shape equality is not strictly
        # necessary as it should be taken care of in np.array_equal, but I'd
        # rather explicitly bail before comparing IDs or data. Use array_equal
        # instead of (a == b).all() because of this issue:
        #     http://stackoverflow.com/a/10582030
        if not isinstance(other, self.__class__):
            eq = False
        elif self.shape != other.shape:
            eq = False
        elif self.sample_ids != other.sample_ids:
            eq = False
        elif not np.array_equal(self.data, other.data):
            eq = False

        return eq

    def to_file(self, out_f, delimiter='\t', memory_efficient=True):
        formatted_sids = self._format_sample_ids(delimiter)
        out_f.write(formatted_sids)
        out_f.write('\n')

        if memory_efficient:
            for sid, vals in izip(self.sample_ids, self.data):
                out_f.write(delimiter.join([sid] + map(str, vals)))
                out_f.write('\n')
        else:
            temp_f = StringIO()
            np.savetxt(temp_f, self.data, delimiter=delimiter, fmt='%s')
            temp_f.seek(0)

            for sid, line in izip(self.sample_ids, temp_f):
                out_f.write(sid)
                out_f.write(delimiter)
                out_f.write(line)
            temp_f.close()

    def _validate(self, data, sample_ids):
        pass

    def _index_list(self, l):
        return dict([(id_, idx) for idx, id_ in enumerate(l)])

    def _format_sample_ids(self, delimiter):
        return delimiter.join([''] + list(self.sample_ids))

    @staticmethod
    def _extract_sample_ids(dm_f, delimiter):
        header_line = None

        for line in dm_f:
            line = line.strip()

            if line and not line.startswith('#'):
                header_line = line
                break

        try:
            dm_f.seek(0)
        except AttributeError:
            pass

        if header_line is None:
            raise MissingHeaderError("Could not find a valid header line "
                    "containing sample IDs in the distance matrix file.")
        else:
            return header_line.split(delimiter)
