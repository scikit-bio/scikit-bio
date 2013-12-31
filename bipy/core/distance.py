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
                          usecols=range(1, len(sample_ids) + 1))
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

    def to_file(self, out_f, delimiter='\t'):
        formatted_sids = self._format_sample_ids(delimiter)

        temp_f = StringIO()
        np.savetxt(temp_f, self.data, delimiter=delimiter, fmt='%s',
                   header=formatted_sids, comments='')

        temp_f.seek(0)
        header = temp_f.readline()
        out_f.write(header)

        for sid, line in izip(self.sample_ids, temp_f):
            out_f.write(sid)
            out_f.write(delimiter)
            out_f.write(line)
        temp_f.close()

    def __str__(self):
        return '%dx%d distance matrix\nSample IDs: %s\n' % (self.shape[0],
                self.shape[1], ', '.join(self.sample_ids)) + str(self.data)

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
