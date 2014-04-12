#!/usr/bin/env python
r"""
Workflow (:mod:`skbio.core.workflow`)
=====================================

.. currentmodule:: skbio.core.workflow

Classes
-------

.. autosummary::
    :toctree: generated/

    Workflow

Examples
--------
>>> from skbio.core.workflow import Workflow
>>> nuc_pattern = 'AATTG'
>>> has_nuc_pattern = lambda s: s[:len(nuc_pattern)] == nuc_pattern
>>> class SequenceProcessor(Workflow):
...    def initialize_state(self, item):
...        # Setup the state for a new item (e.g., a new sequence)
...        self.state = item
...    @method(priority=100)
...    def check_length(self):
...        # Always make sure the sequence is at least 10 nucleotides
...        if len(self.state) < 10:
...            self.failed = True
...    @method(priority=90)
...    @requires(state=has_nuc_pattern)
...    def truncate(self):
...        # Truncate if a specific starting nucleotide pattern is observed
...        self.state = self.state[len(nuc_pattern):]
...    @method(priority=80)
...    @requires(option='reverse', values=True)
...    def reverse(self):
...        # Reverse the sequence if indicatd at runtime
...        self.state = self.state[::-1]
>>> wf = SequenceProcessor(state=None, options={'reverse=': False})
>>> seqs = ['AAAAAAATTTTTTT', 'ATAGACC', 'AATTGCCGGAC', 'ATATGAACAAA']
>>> def success_f(obj):
...     return "SUCCESS: %s" % obj.state
>>> def fail_f(obj):
...     return "FAIL: %s" % obj.state
>>> for result in wf(seqs, success_callback=success_f, fail_callback=fail_f):
...     print result
SUCCESS: AAAAAAATTTTTTT
FAIL: ATAGACC
SUCCESS: CCGGAC
SUCCESS: ATATGAACAAA
>>> wf = SequenceProcessor(state=None, options={'reverse':True}, debug=True)
>>> gen = wf(seqs, fail_callback=lambda x: x.state)
>>> gen.next()
'TTTTTTTAAAAAAA'
>>> print wf.failed
False
>>> gen.next()
'ATAGACC'
>>> print wf.failed
True
>>> print wf.debug_trace
set([('check_length', 0)])
>>> gen.next() #
'CAGGCC'
>>> print wf.failed
False
>>> wf.debug_trace
set([('check_length', 0), ('truncate', 1), ('reverse', 2)])
>>> wf.debug_pre_state[('truncate', 1)]
'AATTGCCGGAC'
>>> wf.debug_post_state[('truncate', 1)]
'CCGGAC'
>>> from skbio.core.workflow import not_none, anything
>>> class Ex(Workflow):
...     @method()
...     @requires(option='foo', values=not_none)
...     def do_something(self):
...         pass
...     @method()
...     @requires(option='bar', values=anything)
...     def do_something_else(self):
...         pass
...     @method()
...     @requires(option='foobar', values=[1,2,3])
...     def do_something_awesome(self):
...         pass
...
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from copy import deepcopy
from time import time
from functools import update_wrapper
from collections import Iterable
from types import MethodType


class NotExecuted(object):
    """Helper object to track if a method was executed"""
    def __init__(self):
        self.msg = None
        self._ghetto_identity = "doc test does insane things"

    def __call__(self, msg):
        self.msg = msg
        return self

    def __eq__(self, other):
        try:
            return self._ghetto_identity == other._ghetto_identity
        except:
            return False

_not_executed = NotExecuted()


class Exists(object):
    """Stub object to assist with ``requires`` when a value exists"""
    def __contains__(self, item):
        return True
anything = Exists()  # external, for when a value can be anything


class NotNone(object):
    def __contains__(self, item):
        if item is None:
            return False
        else:
            return True
not_none = NotNone()


class Workflow(object):
    """Arbitrary workflow support structure

    Methods that are considered to be directly part of the workflow must
    be decorated with ``method``. The workflow methods offer a mechanism to
    logically group functionality together, and are free to make subsequent
    calls to other methods.

    All methods of a subclass of Workflow (those with and without the
    ``method`` decoration) can take advantage of the ``requires`` decorator
    to specify any option or state requirements for the decorated function.

    Parameters
    ----------
    state : object
        State can be anything or nothing. This is dependent on the
        workflow as in some cases, it is useful to preallocate state
        while in other workflows state may be ignored.
    short_circuit : bool
        if True, enables ignoring function methods when a given item
        has failed
    debug : bool
        Enable debug mode
    options : dict
        runtime options, {'option':values}, that the ``requires``
        decorator can interrogate.
    kwargs : dict
        Additional arguments will be added as member variables to self.
        This is handy if additional contextual information is needed by a
        workflow method (e.g., a lookup table).

    Attributes
    ----------
    state
    short_circuit
    debug
    options
    failed

    """

    def __init__(self, state, short_circuit=True, debug=False, options=None,
                 **kwargs):
        r"""Build thy workflow of self"""
        if options is None:
            self.options = {}
        else:
            self.options = options

        self.short_circuit = short_circuit
        self.failed = False
        self.debug = debug
        self.state = state
        self.iter_ = None

        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                raise AttributeError("%s exists in self!" % k)
            setattr(self, k, v)

        if self.debug:
            self._setup_debug()

    def initialize_state(self, item):
        """Initialize state

        This method is called first prior to any other defined workflow method
        with the exception of _setup_debug_trace if self.debug is True

        Parameters
        ----------
        item : anything
            Workflow dependent
        """
        raise NotImplementedError("Must implement this method")

    def _setup_debug(self):
        """Wrap all methods with debug trace support"""
        # ignore all members of the baseclass
        ignore = set(dir(Workflow))

        for attrname in dir(self):
            if attrname in ignore:
                continue

            attr = getattr(self, attrname)

            if isinstance(attr, MethodType):
                setattr(self, attrname, self._debug_trace_wrapper(attr))

    def _all_wf_methods(self):
        """Get all workflow methods

        Methods are sorted by priority
        """
        methods = []
        for item in dir(self):
            obj = getattr(self, item)
            if hasattr(obj, 'priority'):
                methods.append(obj)

        key = lambda x: getattr(x, 'priority')
        methods_sorted = sorted(methods, key=key, reverse=True)

        if self.debug:
            methods_sorted.insert(0, self._setup_debug_trace)

        return methods_sorted

    def _setup_debug_trace(self):
        """Setup a trace

        The trace is per item iterated over by the workflow. Information about
        each method executed is tracked and keyed by:

            (function name, order of execution)

        Order of execution starts from zero. Multiple calls to the same
        function are independent in the trace.

        The following information is tracked:

            debug_trace : set([key])
            debug_runtime : {key: runtime}
            debug_pre_state : {key: deepcopy(Workflow.state)}, state prior to
                method execution
            debug_post_state : {key: deepcopy(Workflow.state)}, state following
                method execution
        """
        self.debug_counter = 0
        self.debug_trace = set()
        self.debug_runtime = {}
        self.debug_pre_state = {}
        self.debug_post_state = {}

    def __call__(self, iter_, success_callback=None, fail_callback=None):
        """Operate on all the data

        This is the processing engine of the workflow. Callbacks are executed
        following applying all workflow methods to an item from ``iter_``
        (unless ``short_cicruit=True`` in which case method execution for an
        item is stopped if ``failed=True``). Callbacks are provided ``self``
        which allows them to examine any aspect of the workflow.

        Parameters
        ----------
        it : an iterator
        success_callback : method to call on a successful item prior to
            yielding. By default, ``self.state`` is yielded.
        fail_callback : method to call on a failed item prior to yielding. By
            default, failures are ignored.

        .. shownumpydoc
        """
        if success_callback is None:
            success_callback = lambda x: x.state

        self.iter_ = iter_
        workflow = self._all_wf_methods()

        for item in self.iter_:
            self.failed = False

            self.initialize_state(item)
            for func in workflow:
                if self.short_circuit and self.failed:
                    break
                else:
                    func()

            if self.failed:
                if fail_callback is not None:
                    yield fail_callback(self)
            else:
                yield success_callback(self)

        self.iter_ = None

    def _debug_trace_wrapper(self, func):
        """Trace a function call"""
        def wrapped():
            """Track debug information about a method execution"""
            if not hasattr(self, 'debug_trace'):
                cls = self.__class__
                raise AttributeError("%s doesn't have debug_trace!" % cls)

            exec_order = self.debug_counter
            name = func.__name__
            key = (name, exec_order)
            pre_state = deepcopy(self.state)

            self.debug_trace.add(key)
            self.debug_counter += 1

            start_time = time()
            if func() == _not_executed:
                self.debug_trace.remove(key)
            else:
                self.debug_runtime[key] = time() - start_time
                self.debug_pre_state[key] = pre_state
                self.debug_post_state[key] = deepcopy(self.state)

        return update_wrapper(wrapped, func)


class method(object):
    """Decorate a function to indicate it is a workflow method

    Parameters
    ----------
    priority : int
        Specify a priority for the method, the higher the value the higher
        the priority. Priorities are relative to a given workflow

    Returns
    -------
    Method
        A decorated method
    """
    highest_priority = sys.maxsize

    def __init__(self, priority=0):
        self.priority = priority

    def __call__(self, func):
        func.priority = self.priority
        return func


class requires(object):
    """Decorator that executes a function if requirements are met

    Parameters
    ----------
    option : any Hashable object
        An option that is required for the decorated method to execute.
        This option will be looked up within the containing ``Workflow``s'
        ``options``,
    values : literally anything
        A required value. This defaults to ``anything`` indicating that
        the only requirement is that the ``option`` exists. It can be
        useful to specify ``not_none`` which indicates that the
        requirement is satisfied if the ``option`` exists and it holds
        a value that is not ``None``. Values also supports iterables
        or singular values
    state : Function
        A requirement on workflow state. This must be a function that
        accepts a single argument, and returns ``True`` to indicate
        the requirement is satisfied, or ``False`` to indicate the
        requirement is not satisfied. This method will be passed the
        containing ``Workflow``s' ``state`` member variable.
    """
    def __init__(self, option=None, values=anything, state=None):
        # self here is the requires object
        self.option = option
        self.required_state = state

        if values is anything:
            self.values = anything
        elif values is not_none:
            self.values = not_none
        elif isinstance(values, set):
            self.values = values
        else:
            if isinstance(values, str):
                self.values = values
            elif isinstance(values, Iterable):
                self.values = set(values)
            else:
                self.values = set([values])

    def __call__(self, func):
        """Wrap a function

        func : the function to wrap
        """
        def decorated(dec_self):
            """A decorated function that has requirements

            dec_self : this is "self" for the decorated function
            """
            if self.required_state is not None:
                if not self.required_state(dec_self.state):
                    return _not_executed

            s_opt = self.option
            ds_opts = dec_self.options

            # if this is a function that doesn't have an option to validate
            if s_opt is None:
                func(dec_self)

            # if the option exists in the Workflow
            elif s_opt in ds_opts:
                val = ds_opts[s_opt]

                # if the value just needs to be not None
                if self.values is not_none and val is not None:
                    func(dec_self)

                # otherwise make sure the value is acceptable
                elif val in self.values:
                    func(dec_self)

                else:
                    return _not_executed

            else:
                return _not_executed

        return update_wrapper(decorated, func)
