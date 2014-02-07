#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from itertools import izip
from bipy.core.workflow import (Workflow, requires, priority,
                                no_requirements, not_none)
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
        return izip(*to_gen)


class MockWorkflow(Workflow):
    def _sanity_check(self):
        pass

    def _allocate_final_state(self):
        self.FinalState = None

    @priority(90)
    @requires(Option='A', Values=True)
    def wf_groupA(self, item):
        self.methodA1(item)
        self.methodA2(item)

    @requires(Option='B', Values=True)
    def wf_groupB(self, item):
        self.methodB1(item)
        self.methodB2(item)

    @priority(10)
    @requires(Option='C', Values=True)
    def wf_groupC(self, item):
        self.methodC1(item)
        self.methodC2(item)

    @requires(IsValid=False)  # always execute
    def methodA1(self, item):
        name = 'A1'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)

    def methodA2(self, item):
        name = 'A2'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)

    @requires(IsValid=False)
    def methodB1(self, item):
        name = 'B1'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
            self.FinalState = 'failed'
        else:
            self.FinalState = (name, item)

    @requires(Option='foo', Values=[1, 2, 3])
    def methodB2(self, item):
        name = 'B2'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
            self.FinalState = 'failed'
        else:
            self.FinalState = (name, item)

    @no_requirements
    def methodC1(self, item):
        name = 'C1'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)

    @requires(IsValid=True, Option='C2', Values=[1, 2, 3])
    def methodC2(self, item):
        name = 'C2'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)


class WorkflowTests(TestCase):
    def setUp(self):
        opts = {'A': True, 'C': True}
        self.obj_short = MockWorkflow(options=opts)
        self.obj_debug = MockWorkflow(debug=True, options=opts)
        self.obj_noshort = MockWorkflow(short_circuit=False, options=opts)

    def test_untagged_wf_method(self):
        class WFTest(Workflow):
            @no_requirements
            def wf_1(self):
                pass

            def wf_2(self):
                pass

        with self.assertRaises(AttributeError):
            _ = WFTest()

    def test_get_workflow_debug(self):
        gen = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})
        exp_wf = [self.obj_debug.wf_setup_debug_trace,
                  self.obj_debug.wf_groupA,
                  self.obj_debug.wf_groupC]
        obs_gen, obs_wf = self.obj_debug._get_workflow(gen)

        self.assertEqual(obs_wf, exp_wf)
        self.assertEqual(list(obs_gen), [1, 2, 3, 4, 5])

        self.assertEqual(self.obj_debug.stats, {})
        self.assertTrue(self.obj_debug.short_circuit)

    def test_debug_trace(self):
        gen = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})
        obj = self.obj_debug(gen)

        exp = ('C1', 1)
        obs = obj.next()
        self.assertEqual(obs, exp)

        exp = ['wf_groupA',
               'methodA1',
               'methodA2',
               'wf_groupC',
               'methodC1',
               'methodC2']
        obs = self.obj_debug.debug_trace
        self.assertEqual(obs, exp)

    def test_get_workflow(self):
        gen = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})
        exp_wf = [self.obj_short.wf_groupA, self.obj_short.wf_groupC]
        obs_gen, obs_wf = self.obj_short._get_workflow(gen)

        self.assertEqual(obs_wf, exp_wf)
        self.assertEqual(list(obs_gen), [1, 2, 3, 4, 5])

        self.assertEqual(self.obj_short.stats, {})
        self.assertTrue(self.obj_short.short_circuit)

    def test_init(self):
        self.assertEqual(self.obj_short.options, {'A': True, 'C': True})
        self.assertEqual(self.obj_short.stats, {})
        self.assertTrue(self.obj_short.short_circuit)
        self.assertEqual(self.obj_noshort.options, {'A': True, 'C': True})
        self.assertEqual(self.obj_noshort.stats, {})
        self.assertFalse(self.obj_noshort.short_circuit)

    def test_all_wf_methods(self):
        # note on priority: groupA:90, groupC:10, groupB:0 (default)
        exp = [self.obj_short.wf_groupA, self.obj_short.wf_groupC,
               self.obj_short.wf_groupB]
        obs = self.obj_short._all_wf_methods()
        self.assertEqual(obs, exp)

    def test_call_AC_no_fail(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})
        sf = lambda x: x.FinalState  # success function

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 5}
        # C2 isn't executed as its requirements aren't met in the options
        exp_result = [('C1', 1), ('C1', 2), ('C1', 3), ('C1', 4), ('C1', 5)]

        obs_result = list(self.obj_short(iter_, sf, None))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(self.obj_short.stats, exp_stats)

    def test_call_AC_fail(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 'fail A2', 4, 5]})
        sf = lambda x: x.FinalState  # success function
        ff = lambda x: x.FinalState  # failed function

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 4, 'C2': 4}

        self.obj_short.options['C2'] = 1
        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_short(iter_, sf, ff)

        r1 = gen.next()
        self.assertEqual(r1, ('C2', 1))
        self.assertFalse(self.obj_short.failed)

        r2 = gen.next()
        self.assertEqual(r2, ('C2', 2))
        self.assertFalse(self.obj_short.failed)

        r3 = gen.next()
        self.assertEqual(self.obj_short.FinalState, ('A2', 'fail A2'))
        self.assertTrue(self.obj_short.failed)
        self.assertEqual(r3, ('A2', 'fail A2'))

        r4 = gen.next()
        self.assertEqual(r4, ('C2', 4))
        self.assertFalse(self.obj_short.failed)

        r5 = gen.next()
        self.assertEqual(r5, ('C2', 5))
        self.assertFalse(self.obj_short.failed)

        self.assertEqual(self.obj_short.stats, exp_stats)

    def test_call_AC_fail_noshort(self):
        iter_ = construct_iterator(**{'iter_x': [1, 2, 'fail A2', 4, 5]})
        sf = lambda x: x.FinalState  # success function
        ff = lambda x: x.FinalState  # failed function

        exp_stats = {'A1': 5, 'A2': 5, 'C1': 5}

        # pass in a failed callback to capture the result, and pause execution
        gen = self.obj_noshort(iter_, sf, ff)

        r1 = gen.next()
        self.assertEqual(r1, ('C1', 1))
        self.assertFalse(self.obj_noshort.failed)

        r2 = gen.next()
        self.assertEqual(r2, ('C1', 2))
        self.assertFalse(self.obj_noshort.failed)

        _ = gen.next()
        self.assertEqual(self.obj_noshort.FinalState, ('C1', 'fail A2'))
        self.assertTrue(self.obj_noshort.failed)

        r4 = gen.next()
        self.assertEqual(r4, ('C1', 4))
        self.assertFalse(self.obj_noshort.failed)

        r5 = gen.next()
        self.assertEqual(r5, ('C1', 5))
        self.assertFalse(self.obj_noshort.failed)

        self.assertEqual(self.obj_noshort.stats, exp_stats)


class MockWorkflowReqTest(Workflow):
    def _sanity_check(self):
        pass

    def _allocate_final_state(self):
        self.FinalState = None

    @priority(5)
    @requires(ValidData=lambda x: x < 3)
    def wf_needs_data(self, item):
        name = 'needs_data'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)

    @priority(10)
    @no_requirements
    def wf_always_run(self, item):
        name = 'always_run'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)

    @priority(20)
    @requires(Option='cannot_be_none', Values=not_none)
    def wf_run_if_not_none(self, item):
        name = 'run_if_not_none'
        self.stats[name] += 1
        if item == 'fail %s' % name:
            self.failed = True
        self.FinalState = (name, item)


class RequiresTests(TestCase):
    def test_validdata(self):
        obj = MockWorkflowReqTest()
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5}
        exp_result = [('needs_data', 1), ('needs_data', 2), ('always_run', 3),
                      ('always_run', 4), ('always_run', 5)]

        obs_result = list(obj(single_iter))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_not_none_avoid(self):
        obj = MockWorkflowReqTest({'cannot_be_none': None})
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5}
        exp_result = [('needs_data', 1), ('needs_data', 2), ('always_run', 3),
                      ('always_run', 4), ('always_run', 5)]

        obs_result = list(obj(single_iter))

        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_not_none_execute(self):
        obj = MockWorkflowReqTest(options={'cannot_be_none': True}, debug=True)
        single_iter = construct_iterator(**{'iter_x': [1, 2, 3, 4, 5]})

        exp_stats = {'needs_data': 2, 'always_run': 5, 'run_if_not_none': 5}
        exp_result = [('needs_data', 1), ('needs_data', 2), ('always_run', 3),
                      ('always_run', 4), ('always_run', 5)]

        obs_result = list(obj(single_iter))
        self.assertEqual(obs_result, exp_result)
        self.assertEqual(obj.stats, exp_stats)

    def test_methodb1(self):
        obj = MockWorkflow()
        obj.methodB1('test')
        self.assertEqual(obj.FinalState, ('B1', 'test'))
        self.assertFalse(obj.failed)

        # methodb1 executes regardless of if self.failed
        obj.failed = True
        obj.methodB1('test 2')
        self.assertEqual(obj.FinalState, ('B1', 'test 2'))

        obj.failed = False
        obj.methodB1('fail B1')
        self.assertEqual(obj.FinalState, 'failed')

        self.assertEqual(obj.stats, {'B1': 3})

    def test_methodb2_accept(self):
        # methodb2 is setup to be valid when foo is in [1,2,3], make sure we
        # can execute
        obj = MockWorkflow(options={'foo': 1})
        obj.methodB2('test')
        self.assertEqual(obj.FinalState, ('B2', 'test'))
        self.assertEqual(obj.stats, {'B2': 1})

        # methodb2 will not execute if self.failed
        obj.failed = True
        obj.methodB2('test 2')
        self.assertEqual(obj.FinalState, ('B2', 'test'))
        self.assertEqual(obj.stats, {'B2': 1})

    def test_methodb2_ignore(self):
        # methodb2 is setup to be valid when foo is in [1, 2, 3], make sure
        # we do not execute
        obj = MockWorkflow(options={'foo': 'bar'})
        obj.methodB2('test')
        self.assertEqual(obj.FinalState, None)
        self.assertEqual(obj.stats, {})


class PriorityTests(TestCase):
    def test_dec(self):
        @priority(10)
        def foo(x, y, z):
            """doc check"""
            return x + y + z

        self.assertEqual(foo.Priority, 10)
        self.assertEqual(foo.__name__, 'foo')
        self.assertEqual(foo.__doc__, 'doc check')

if __name__ == '__main__':
    main()
