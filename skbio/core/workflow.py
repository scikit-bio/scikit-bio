#!/usr/bin/env python

"""Perform multiple method calls, determined at runtime, on independent items

Construct arbitrarily complex workflows in which the specific methods run are
determined at runtime. These methods are applied to items that are assumed to
be independent.

As an example:

class MyWorkflow(Workflow):
    def initialize_state(self, item):
        self.state = item

    @Workflow.method(priority=100)
    def wf_mul(self):
        self.state *= self.state

    @Workflow.method(priority=10)
    @Workflow.requires(option='double')
    def wf_double(self):
        self.state += self.state

    @Workflow.method()
    @Workflow.requires(option='sub_value', values=[1,5,10])
    def wf_sub(self):
        self.state -= self.options['sub_value']


# ((i * i) * 2) - 5
wf = MyWorkflow(0, options={'double':None, 'sub_value':5})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i) - 10
wf = MyWorkflow(0, options={'sub_value':10})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i)
wf = MyWorkflow(0)
gen = (i for i in range(10))
for i in wf(gen):
    print i
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

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

    def __call__(self, msg):
        self.msg = msg
        return self
_not_executed = NotExecuted()


class Exists(object):
    """Stub object to assist with Workflow.requires when a value exists"""
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
    """Arbitrary workflow support structure"""

    def __init__(self, state, short_circuit=True, debug=False, options=None,
                 **kwargs):
        """Build thy self

        state : an allocated object to store state
        short_circuit : if True, enables ignoring function methods when a given
            item has failed
        debug : Enable debug mode
        options : runtime options, {'option':values}
        kwargs : Additional arguments will be added to self

        Methods that are considered to be directly part of the workflow must
        be decorated with Workflow.method. The workflow methods offer a
        mechanism to logically group functionality together, and are free to
        make subsequent calls to other methods.

        All methods of a subclass of Workflow (those with and without the
        Workflow.method decoration) can take advantage of the Workflow.requires
        decorator to specify any option or state requirements for the
        decorated function.
        """
        if options is None:
            self.options = {}
        else:
            self.options = options

        self.short_circuit = short_circuit
        self.failed = False
        self.debug = debug
        self.state = state

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

        it : an iterator
        success_callback : method to call on a successful item prior to
            yielding
        fail_callback : method to call on a failed item prior to yielding
        """
        if success_callback is None:
            success_callback = lambda x: x.state

        workflow = self._all_wf_methods()

        for item in iter_:
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

    ### Decorators ###

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
            if func() is _not_executed:
                self.debug_trace.remove(key)
            else:
                self.debug_runtime[key] = time() - start_time
                self.debug_pre_state[key] = pre_state
                self.debug_post_state[key] = deepcopy(self.state)

        return update_wrapper(wrapped, func)

    class method(object):
        """Decorate a function to indicate it is a workflow method"""
        highest_priority = sys.maxsize

        def __init__(self, priority=0):
            self.priority = priority

        def __call__(self, func):
            func.priority = self.priority
            return func

    class requires(object):
        """Decorator that executes a function if requirements are met"""
        def __init__(self, option=None, values=anything, state=None):
            """
            option : a required option
            values : required values associated with an option
            state : If state is not None and is a function, that function will
                check if state(Workflow.state) is True, or valid, for each
                item in the iterable operated on by Workflow.__call__.
            """
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
