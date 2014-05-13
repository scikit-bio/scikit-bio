#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio.math.stats.spatial import (procrustes, _get_disparity, _center,
                                      _normalize)


class ProcrustesTests(TestCase):

    """test the procrustes module, using floating point numpy arrays
    """

    def setUp(self):
        """creates inputs"""
        # an L
        self.data1 = np.array([[1, 3], [1, 2], [1, 1], [2, 1]], 'd')

        # a larger, shifted, mirrored L
        self.data2 = np.array([[4, -2], [4, -4], [4, -6], [2, -6]], 'd')

        # an L shifted up 1, right 1, and with point 4 shifted an extra .5
        # to the right
        # pointwise distance disparity with data1: 3*(2) + (1 + 1.5^2)
        self.data3 = np.array([[2, 4], [2, 3], [2, 2], [3, 2.5]], 'd')

        # data4, data5 are standardized (trace(A*A') = 1).
        # procrustes should return an identical copy if they are used
        # as the first matrix argument.
        shiftangle = np.pi / 8
        self.data4 = np.array([[1, 0], [0, 1], [-1, 0],
                              [0, -1]], 'd') / np.sqrt(4)
        self.data5 = np.array([[np.cos(shiftangle), np.sin(shiftangle)],
                              [np.cos(np.pi / 2 - shiftangle),
                               np.sin(np.pi / 2 - shiftangle)],
                              [-np.cos(shiftangle),
                               -np.sin(shiftangle)],
                              [-np.cos(np.pi / 2 - shiftangle),
                               -np.sin(np.pi / 2 - shiftangle)]],
                              'd') / np.sqrt(4)

    def test_procrustes(self):
        """tests procrustes' ability to match two matrices.

        the second matrix is a rotated, shifted, scaled, and mirrored version
        of the first, in two dimensions only
        """
        # can shift, mirror, and scale an 'L'?
        a, b, disparity = procrustes(self.data1, self.data2)
        np.testing.assert_allclose(b, a)
        np.testing.assert_almost_equal(disparity, 0.)

        # if first mtx is standardized, leaves first mtx unchanged?
        m4, m5, disp45 = procrustes(self.data4, self.data5)
        np.testing.assert_equal(m4, self.data4)

        # at worst, data3 is an 'L' with one point off by .5
        m1, m3, disp13 = procrustes(self.data1, self.data3)
        self.assertTrue(disp13 < 0.5 ** 2)

    def test_procrustes2(self):
        """procrustes disparity should not depend on order of matrices"""
        m1, m3, disp13 = procrustes(self.data1, self.data3)
        m3_2, m1_2, disp31 = procrustes(self.data3, self.data1)
        np.testing.assert_almost_equal(disp13, disp31)

        # try with 3d, 8 pts per
        rand1 = np.array([[2.61955202,  0.30522265,  0.55515826],
                         [0.41124708, -0.03966978, -0.31854548],
                         [0.91910318,  1.39451809, -0.15295084],
                         [2.00452023,  0.50150048,  0.29485268],
                         [0.09453595,  0.67528885,  0.03283872],
                         [0.07015232,  2.18892599, -1.67266852],
                         [0.65029688,  1.60551637,  0.80013549],
                         [-0.6607528,  0.53644208,  0.17033891]])

        rand3 = np.array([[0.0809969,  0.09731461, -0.173442],
                         [-1.84888465, -0.92589646, -1.29335743],
                         [0.67031855, -1.35957463,  0.41938621],
                         [0.73967209, -0.20230757,  0.52418027],
                         [0.17752796,  0.09065607,  0.29827466],
                         [0.47999368, -0.88455717, -0.57547934],
                         [-0.11486344, -0.12608506, -0.3395779],
                         [-0.86106154, -0.28687488,  0.9644429]])
        res1, res3, disp13 = procrustes(rand1, rand3)
        res3_2, res1_2, disp31 = procrustes(rand3, rand1)
        np.testing.assert_almost_equal(disp13, disp31)

    def test_get_disparity(self):
        """tests get_disparity"""
        disp = _get_disparity(self.data1, self.data3)
        disp2 = _get_disparity(self.data3, self.data1)
        np.testing.assert_equal(disp, disp2)
        np.testing.assert_equal(disp, (3. * 2. + (1. + 1.5 ** 2)))

        d1 = np.append(self.data1, self.data1, 0)
        d3 = np.append(self.data3, self.data3, 0)

        disp3 = _get_disparity(d1, d3)
        disp4 = _get_disparity(d3, d1)
        np.testing.assert_equal(disp3, disp4)
        # 2x points in same configuration should give 2x disparity
        np.testing.assert_equal(disp3, 2. * disp)

    def test_center(self):
        centered_mtx = _center(self.data1)
        column_means = centered_mtx.mean(0)
        for col_mean in column_means:
            np.testing.assert_equal(col_mean, 0.)

    def test_normalize(self):
        norm_mtx = _normalize(self.data1)
        np.testing.assert_equal(np.trace(np.dot(norm_mtx,
                                                np.transpose(norm_mtx))), 1.)

    # match_points isn't yet tested, as it's almost a private function
    # and test_procrustes() tests it implicitly.


if __name__ == '__main__':
    main()
