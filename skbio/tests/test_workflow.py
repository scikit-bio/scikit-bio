# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict
from skbio.workflow import (Exists, NotExecuted, NotNone, Workflow, not_none,
                            requires, method)
from unittest import TestCase, main


def construct_iterator(**kwargs):
    """make an iterator for testing purposes"""
    to_gen = []
    for k in sorted(kwargs):
        if k.startswith('iter'):
            to_gen.append(kwargs[k])
    if len(to_gen) == 1:
        return (x for x in to_gen[0])
    else:
        return zip(*to_gen)


class MockWorkflow(Workflow):
    def initialize_state(self, item):
        self.state[0] = None
        self.state[1] = item

    @method(priority=90)
    @requires(option='A', values=True)
    def wf_groupA(self):
        self.methodA1()
        self.methodA2()

    @method()
    @requires(option='B', values=True)
    def wf_groupB(self):
        self.methodB1()
        self.methodB2()

    @method(priority=10)
    @requires(option='C', values=True)
    def wf_groupC(self):
        self.methodC1()
        self.methodC2()

    def methodA1(self):
        name = 'A1'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]

    def methodA2(self):
        name = 'A2'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]

    def methodB1(self):
        name = 'B1'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
            self.state = 'failed'
        else:
            self.state = [name, self.state[-1]]

    @requires(option='foo', values=[1, 2, 3])
    def methodB2(self):
        name = 'B2'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
            self.state = 'failed'
        else:
            self.state = [name, self.state[-1]]

    def methodC1(self):
        name = 'C1'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]

    @requires(option='C2', values=[1, 2, 3])
    def methodC2(self):
        name = 'C2'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]


class WorkflowTests(TestCase):
    def setUp(self):
        opts = {'A': True, 'C': True}
        self.obj_short = MockWorkflow([None, None], options=opts,
                                      stats=defaultdict(int))
        self.obj_debug = MockWorkflow([None, None], debug=True, options=opts,
                                      stats=defaultdict(int))
        self.obj_noshort = MockWorkflow([None, None], short_circuit=False,
                                        options=opts,
                                        stats=defaultdict(int))

    def test_debug_trace(self):
        gen = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})
        obj = self.obj_debug(gen)

        exp = ['C1', 1]
        obs = next(obj)
        self.assertEqual(obs, exp)

        exp_trace = set([('wf_groupA', 0),
                         ('methodA1', 1),
                         ('methodA2', 2),
                         ('wf_groupC', 3),
                         ('methodC1', 4)])

        exp_pre_state = {('wf_groupA', 0): [None, 1],
                         ('methodA1', 1): [None, 1],
                         ('methodA2', 2): ['A1', 1],
                         ('wf_groupC', 3): ['A2', 1],
                         ('methodC1', 4): ['A2', 1]}

        exp_post_state = {('wf_groupA', 0): ['A2', 1],
                          ('methodA1', 1): ['A1', 1],
                          ('methodA2', 2): ['A2', 1],
                          ('wf_groupC', 3): ['C1', 1],
                          ('methodC1', 4): ['C1', 1]}

        obs_trace = self.obj_debug.debug_trace
        obs_pre_state = self.obj_debug.debug_pre_state
        obs_post_state = self.obj_debug.debug_post_state

        self.assertEqual(obs_trace, exp_trace)
        self.assertEqual(obs_pre_state, exp_pre_state)
        self.assertEqual(obs_post_state, exp_post_state)

    def test_init(self):
        self.assertEqual(self.obj_short.options, {'A': True, 'C': True})
        self.assertEqual(self.obj_short.stats, {})
        self.assertTrue(self.obj_short.short_circuit)
        self.assertEqual(self.obj_noshort.options, {'A': True, 'C': True})
        self.assertEqual(self.obj_noshort.stats, {})
        self.assertFalse(self.obj_noshort.short_circuit)

    def test_init_reserved_attributes(self):
        with self.assertRaises(AttributeError):
            Workflow('foo', failed=True)

    def test_all_wf_methods(self):
        # note on priority: groupA:90, groupC:10, groupB:0 (default)
        exp = [self.obj_short.wf_groupA, self.obj_short.wf_groupC,
               self.obj_short.wf_groupB]
        obs = self.obj_short._all_wf_methods()
        self.assertEqual(obs, exp)

    def test_call_AC_no_fail(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        # success function
        def sf(x):
            return x.state[:]

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 5}
        # C2 isn't executed as its requirements aren't met in the options
        exp_result = [['C1', 1], ['C1', 2], ['C1', 3], ['C1', 4], ['C1', 5]]

        obs_result = list(self.obj_short(iter_, sf, None))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(self.obj_short.stats, exp_stats)

    def test_call_AC_fail(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 'fail A2', 4, 5]})

        # success function
        def sf(x):
            return x.state[:]

        ff = sf  # failed function

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 4, 'C2': 4}

        self.obj_short.options['C2'] = 1
        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_short(iter_, sf, ff)

        r1 = next(gen)
        self.assertEqual(r1, ['C2', 1])
        self.assertFalse(self.obj_short.failed)

        r2 = next(gen)
        self.assertEqual(r2, ['C2', 2])
        self.assertFalse(self.obj_short.failed)

        r3 = next(gen)
        self.assertEqual(self.obj_short.state, ['A2', 'fail A2'])
        self.assertTrue(self.obj_short.failed)
        self.assertEqual(r3, ['A2', 'fail A2'])

        r4 = next(gen)
        self.assertEqual(r4, ['C2', 4])
        self.assertFalse(self.obj_short.failed)

        r5 = next(gen)
        self.assertEqual(r5, ['C2', 5])
        self.assertFalse(self.obj_short.failed)

        self.assertEqual(self.obj_short.stats, exp_stats)

    def test_call_AC_fail_noshort(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 'fail A2', 4, 5]})

        # success function
        def sf(x):
            return x.state[:]

        ff = sf  # failed function

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 5}

        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_noshort(iter_, sf, ff)

        r1 = next(gen)
        self.assertEqual(r1, ['C1', 1])
        self.assertFalse(self.obj_noshort.failed)

        r2 = next(gen)
        self.assertEqual(r2, ['C1', 2])
        self.assertFalse(self.obj_noshort.failed)

        next(gen)
        self.assertEqual(self.obj_noshort.state, ['C1', 'fail A2'])
        self.assertTrue(self.obj_noshort.failed)

        r4 = next(gen)
        self.assertEqual(r4, ['C1', 4])
        self.assertFalse(self.obj_noshort.failed)

        r5 = next(gen)
        self.assertEqual(r5, ['C1', 5])
        self.assertFalse(self.obj_noshort.failed)

        self.assertEqual(self.obj_noshort.stats, exp_stats)


class MockWorkflowReqTest(Workflow):
    def _allocate_state(self):
        self.state = None

    def initialize_state(self, item):
        self.state = [None, item]

    @method(priority=5)
    @requires(state=lambda x: x[-1] < 3)
    def wf_needs_data(self):
        name = 'needs_data'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]

    @method(priority=10)
    def wf_always_run(self):
        name = 'always_run'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]

    @method(priority=20)
    @requires(option='cannot_be_none', values=not_none)
    def wf_run_if_not_none(self):
        name = 'run_if_not_none'
        self.stats[name] += 1
        if self.state[-1] == 'fail %s' % name:
            self.failed = True
        self.state = [name, self.state[-1]]


class RequiresTests(TestCase):
    def test_validdata(self):
        obj = MockWorkflowReqTest([None, None], stats=defaultdict(int))
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5}
        exp_result = [['needs_data', 1], ['needs_data', 2], ['always_run', 3],
                      ['always_run', 4], ['always_run', 5]]

        obs_result = list(obj(single_iter))
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_not_none_avoid(self):
        obj = MockWorkflowReqTest([None, None], {'cannot_be_none': None},
                                  stats=defaultdict(int))
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5}
        exp_result = [['needs_data', 1], ['needs_data', 2], ['always_run', 3],
                      ['always_run', 4], ['always_run', 5]]

        obs_result = list(obj(single_iter))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_not_none_execute(self):
        obj = MockWorkflowReqTest([None, None],
                                  options={'cannot_be_none': True}, debug=True,
                                  stats=defaultdict(int))
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5, 'run_if_not_none': 5}
        exp_result = [['needs_data', 1], ['needs_data', 2], ['always_run', 3],
                      ['always_run', 4], ['always_run', 5]]

        obs_result = list(obj(single_iter))
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_methodb1(self):
        obj = MockWorkflow([None, None], stats=defaultdict(int))
        obj.initialize_state('test')
        obj.methodB1()
        self.assertEqual(obj.state, ['B1', 'test'])
        self.assertFalse(obj.failed)

        # methodb1 executes regardless of if self.failed
        obj.failed = True
        obj.initialize_state('test 2')
        obj.methodB1()
        self.assertEqual(obj.state, ['B1', 'test 2'])

        obj.failed = False
        obj.state = [None, 'fail B1']
        obj.methodB1()
        self.assertEqual(obj.state, 'failed')

        self.assertEqual(obj.stats, {'B1': 3})

    def test_methodb2_accept(self):
        # methodb2 is setup to be valid when foo is in [1,2,3], make sure we
        # can execute
        obj = MockWorkflow([None, None], options={'foo': 1},
                           stats=defaultdict(int))
        obj.initialize_state('test')
        obj.methodB2()
        self.assertEqual(obj.state, ['B2', 'test'])
        self.assertEqual(obj.stats, {'B2': 1})

    def test_methodb2_ignore(self):
        # methodb2 is setup to be valid when foo is in [1, 2, 3], make sure
        # we do not execute
        obj = MockWorkflow([None, None], options={'foo': 'bar'},
                           stats=defaultdict(int))
        obj.methodB2()
        self.assertEqual(obj.state, [None, None])
        self.assertEqual(obj.stats, {})


class PriorityTests(TestCase):
    def test_dec(self):
        @method(priority=10)
        def foo(x, y, z):
            """doc check"""
            return x + y + z

        self.assertEqual(foo.priority, 10)
        self.assertEqual(foo.__name__, 'foo')
        self.assertEqual(foo.__doc__, 'doc check')


class NotExecutedTests(TestCase):
    def test_call(self):
        ne = NotExecuted()
        obs = ne('foo')
        self.assertTrue(obs is ne)
        self.assertEqual(obs.msg, 'foo')


class ExistsTests(TestCase):
    def test_contains(self):
        e = Exists()
        self.assertTrue('foo' in e)
        self.assertTrue(None in e)


class NotNoneTests(TestCase):
    def test_contains(self):
        nn = NotNone()
        self.assertTrue('foo' in nn)
        self.assertFalse(None in nn)


if __name__ == '__main__':
    main()
