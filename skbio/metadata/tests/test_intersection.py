# checklist.py:CopyrightHeadersValidator IGNORE
# ----------------------------------------------------------------------------
# This code is taken from bx-python project and added with a test for a new
# function `update` from line 195 to 211. The license for this code is
# included in licenses/bx_python.txt.
# ----------------------------------------------------------------------------

import sys
import unittest
import random

from skbio.metadata._intersection import IntervalObj
from skbio.metadata._intersection import IntervalNode
from skbio.metadata._intersection import IntervalTree


# Note: More bounds checking on input parameters are performed
# within test_interval.py
class NeighborTestCase(unittest.TestCase):

    def setUp(self):
        iv = IntervalNode(50, 59, IntervalObj(50, 59))
        for i in range(0, 110, 10):
            if i == 50:
                continue
            f = IntervalObj(i, i + 9)
            iv = iv.insert(f.start, f.end, f)
        self.intervals = iv

    def test_left(self):
        iv = self.intervals
        self.assertEqual(str(iv.left(60, n=2)),
                         str([IntervalObj(50, 59), IntervalObj(40, 49)]))

        for i in range(10, 100, 10):
            r = iv.left(i, max_dist=10, n=1)
            self.assertEqual(r[0].end,  i - 1)

    def test_toomany(self):
        iv = self.intervals
        self.assertEqual(len(iv.left(60, n=200)), 6)

    def test_right(self):
        iv = self.intervals
        self.assertEqual(str(iv.left(60, n=2)),
                         str([IntervalObj(50, 59), IntervalObj(40, 49)]))

        def get_right_start(b10):
            r = iv.right(b10+1, n=1)
            assert len(r) == 1
            return r[0].start

        for i in range(10, 100, 10):
            self.assertEqual(get_right_start(i), i + 10)

        for i in range(0, 100, 10):
            r = iv.right(i-1, max_dist=10, n=1)
            self.assertEqual(r[0].start, i)


class UpDownStreamTestCase(unittest.TestCase):

    def setUp(self):
        iv = IntervalTree()
        iv.add_interval(IntervalObj(50, 59))
        for i in range(0, 110, 10):
            if i == 50:
                continue
            f = IntervalObj(i, i + 9)
            iv.add_interval(f)
        self.intervals = iv

    def test_upstream(self):
        iv = self.intervals
        upstreams = iv.upstream_of_interval(IntervalObj(59, 60),
                                            num_intervals=200)

        for u in upstreams:
            self.assertTrue(u.end < 59)

        upstreams = iv.upstream_of_interval(IntervalObj(60, 70, strand=-1),
                                            num_intervals=200)
        for u in upstreams:
            self.assertTrue(u.start > 70)

        upstreams = iv.upstream_of_interval(IntervalObj(58, 58, strand=-1),
                                            num_intervals=200)
        for u in upstreams:
            self.assertTrue(u.start > 59)

    def test_downstream(self):
        iv = self.intervals
        downstreams = iv.downstream_of_interval(IntervalObj(59, 60),
                                                num_intervals=200)
        for d in downstreams:
            self.assertTrue(d.start > 60)

        downstreams = iv.downstream_of_interval(IntervalObj(59, 60, strand=-1),
                                                num_intervals=200)
        for d in downstreams:
            self.assertTrue(d.start < 59)

    def test_n(self):
        iv = self.intervals
        for i in range(0, 90, 10):
            r = iv.after(i, max_dist=20, num_intervals=2)
            self.assertEqual(r[0].start, i + 10)
            self.assertEqual(r[1].start, i + 20)

            r = iv.after_interval(IntervalObj(i, i),
                                  max_dist=20, num_intervals=2)
            self.assertEqual(r[0].start, i + 10)
            self.assertEqual(r[1].start, i + 20)


class LotsaTestCase(unittest.TestCase):
    """ put lotsa data in the tree and make sure it works"""
    def setUp(self):
        iv = IntervalNode(1, 2, IntervalObj(1, 2))
        self.max = 1000000
        for i in range(0, self.max, 10):
            f = IntervalObj(i, i)
            iv = iv.insert(f.start, f.end, f)

        for i in range(600):
            iv = iv.insert(0, 1, IntervalObj(0, 1))
        self.intervals = iv

    def test_count(self):
        iv = self.intervals

        r = iv.right(1, n=33)
        self.assertEqual(len(r), 33)

        left = iv.left(1, n=33)
        self.assertEqual(len(left), 1)

        u = iv.right(1, n=9999)
        self.assertEqual(len(u), 250)

        # now increase max_dist
        u = iv.right(1, n=9999, max_dist=99999)
        self.assertEqual(len(u), 9999)

    def test_max_dist(self):
        iv = self.intervals
        r = iv.right(1, max_dist=0, n=10)
        self.assertEqual(len(r), 0)

        for n, d in enumerate(range(10, 1000, 10)):
            r = iv.right(1, max_dist=d, n=10000)
            self.assertEqual(len(r), n + 1)

    def test_find(self):
        iv = self.intervals
        path = sys.path[:]
        sys.path = sys.path[2:]
        sys.path = path
        for t in range(25):
            start = random.randint(0, self.max - 10000)
            end = start + random.randint(100, 10000)
            results = iv.find(start, end)
            for feat in results:
                self.assertTrue(
                        (feat.end >= start and feat.end <= end) or
                        (feat.start <= end and feat.start >= start))


class IntervalTreeTest(unittest.TestCase):
    def setUp(self):
        iv = IntervalTree()
        n = 0
        for i in range(1, 1000, 80):
            iv.insert(i, i + 10, dict(value=i*i))
            # add is synonym for insert.
            iv.add(i + 20, i + 30, dict(astr=str(i*i)))

            # or insert/add an interval object with start, end attrs.
            iv.insert_interval(IntervalObj(i + 40, i + 50,
                                           value=dict(astr=str(i*i))))
            iv.add_interval(IntervalObj(i + 60, i + 70,
                                        value=dict(astr=str(i*i))))

            n += 4
        self.intervals = self.iv = iv
        self.nintervals = n

    def test_find(self):
        r = self.iv.find(100, 200)
        self.assertEqual(len(r), 5)

    def test_edge_cases(self):
        iv = IntervalTree()
        iv.insert(1, 1, 'foo')
        iv.insert(3, 7, 'spam')
        iv.insert(8, 8, 'abc')
        self.assertEqual(iv.find(0, 1), [])
        self.assertEqual(iv.find(1, 1), ['foo'])
        self.assertEqual(iv.find(1, 2), ['foo'])
        self.assertEqual(iv.find(2, 3), [])
        self.assertEqual(iv.find(3, 3), ['spam'])
        self.assertEqual(iv.find(3, 4), ['spam'])
        self.assertEqual(iv.find(6, 7), ['spam'])
        self.assertEqual(iv.find(7, 7), [])
        self.assertEqual(iv.find(0, 8), ['foo', 'spam'])
        self.assertEqual(iv.find(8, 8), ['abc'])
        self.assertEqual(iv.find(8, 9), ['abc'])
        self.assertEqual(iv.find(9, 9), [])
        self.assertEqual(iv.find(0, 10), ['foo', 'spam', 'abc'])
        self.assertEqual(iv.find(6, 10), ['spam', 'abc'])
        self.assertEqual(iv.find(7, 9), ['abc'])

    def test_traverse(self):
        a = []
        fn = a.append

        self.iv.traverse(fn)
        self.assertEqual(len(a), self.nintervals)

    def test_empty(self):
        iv = IntervalTree()
        self.assertEqual([], iv.find(100, 300))
        self.assertEqual([], iv.after(100))
        self.assertEqual([], iv.before(100))
        self.assertEqual([], iv.after_interval(100))
        self.assertEqual([], iv.before_interval(100))
        self.assertEqual([], iv.upstream_of_interval(100))
        self.assertEqual([], iv.downstream_of_interval(100))
        self.assertEqual(None, iv.traverse(lambda x: x.append(1)))

    def test_public_interval(self):
        self.iv.traverse(lambda ival: self.assertTrue(ival.interval))

    def test_update(self):
        i = 1
        self.iv.update(i, i + 10, dict(value=i*i), dict(value=-1))
        self.assertEqual([dict(value=-1)], self.iv.find(i, i+10))


if __name__ == "__main__":
    unittest.main()
