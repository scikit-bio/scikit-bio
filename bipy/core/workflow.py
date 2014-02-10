#!/usr/bin/env python

"""Perform multiple method calls, determined at runtime, on independent items

Construct arbitrarily complex workflows in which the specific methods run are
determined at runtime. These methods are applied to items that are assumed to
be independent.

As an example:

class MyWorkflow(Workflow):
    def _allocate_final_state(self):
        self.final_state = None

    def _sanity_check(self):
        pass

    @priority(100)
    @no_requirements
    def wf_mul(self, item):
        self.final_state *= item

    @priority(10)
    @requires(option='add_value')
    def wf_add(self, item):
        self.final_state += item

    @requires(option='sub_value', values=[1,5,10])
    def wf_sub(self, item):
        self.final_state -= item
        self.final_state -= self.options['sub_value']

    @priority(1000)
    @requires(valid_state=False)
    def wf_init(self, item):
        self.final_state = item

# (i * i) + i - i - 5
wf = MyWorkflow(options={'add_value':None, 'sub_value':5})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i) - i - 10
wf = MyWorkflow(options={'sub_value':10})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i)
wf = MyWorkflow()
gen = (i for i in range(10))
for i in wf(gen):
    print i
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys
from itertools import chain
from functools import update_wrapper
from collections import Iterable, defaultdict
from types import MethodType

# thank you Flask project...
_executed = object()  # internal, tag for an executed method
not_none = object()   # external, for when a value can be anything except None


class Exists(object):
    def __contains__(self, item):
        return True
anything = Exists()  # external, for when a value can be anything


def _debug_trace_wrapper(obj, f):
    """Trace a function call"""
    def wrapped(self, *args, **kwargs):
        if not hasattr(obj, 'debug_trace'):
            cls = obj.__class__
            raise AttributeError("%s doesn't have debug_trace!" % cls)

        obj.debug_trace.append(f.__name__)
        return f(self, *args, **kwargs)

    return update_wrapper(wrapped, f)


def _tag_function(f):
    """Tag, you're it"""
    setattr(f, '_workflow_tag', None)


class priority(object):
    """Decorate a function priority"""
    highest = sys.maxint

    def __init__(self, priority):
        self.priority = priority

    def __call__(self, f):
        f.priority = self.priority
        return f


def no_requirements(f):
    """Decorate a function to indicate there are no requirements"""
    def decorated(self, *args, **kwargs):
        """Simply execute the function"""
        f(self, *args, **kwargs)
        return _executed
    _tag_function(decorated)
    return update_wrapper(decorated, f)


class requires(object):
    """Decorator that executes a function if requirements are met"""
    def __init__(self, valid_state=True, option=None, values=anything,
                 valid_data=None):
        """
        valid_state : execute the function if self.failed is False
        option : a required option
        values : required values associated with an option
        valid_data : data level requirements, this must be a function with the
            following signature: f(*args, **kwargs) returning True. NOTE: if
            valid_data returns False on the first item evaluated, the decorated
            function will be removed from the remaining workflow.
        """
        # self here is the requires object
        self.valid_state = valid_state
        self.option = option
        self.valid_data = valid_data

        if values is anything:
            self.values = anything
        elif values is not_none:
            self.values = not_none
        elif not isinstance(values, set):
            if isinstance(values, str):
                self.values = values
            elif isinstance(values, Iterable):
                self.values = set(values)
            else:
                self.values = set([values])
        else:
            self.values = values

    def do_short_circuit(self, wrapped):
        return self.valid_state and (wrapped.failed and wrapped.short_circuit)

    def __call__(self, f):
        """Wrap a function

        f : the function to wrap
        """
        ### not sure how I feel about having multiple functions in here.
        ### also, the handling of Data is a bit dirty as it is now replicated
        ### over these functions. It is ideal to keep the functions slim, thus
        ### the multiple functions, but this could explode if not careful
        def decorated_with_option(dec_self, *args, **kwargs):
            """A decorated function that has an option to validate

            dec_self : this is "self" for the decorated function
            """
            if self.do_short_circuit(dec_self):
                return

            if self.valid_data is not None:
                if not self.valid_data(*args, **kwargs):
                    return

            s_opt = self.option
            ds_opts = dec_self.options

            # if the option exists in the Workflow
            if s_opt in ds_opts:
                v = ds_opts[s_opt]

                # if the value just needs to be not None
                if self.values is not_none and v is not None:
                    f(dec_self, *args, **kwargs)
                    return _executed

                # otherwise make sure the value is acceptable
                elif v in self.values:
                    f(dec_self, *args, **kwargs)
                    return _executed

        def decorated_without_option(dec_self, *args, **kwargs):
            """A decorated function that does not have an option to validate

            dec_self : this is "self" for the decorated function
            """
            if self.do_short_circuit(dec_self):
                return

            if self.valid_data is not None:
                if not self.valid_data(*args, **kwargs):
                    return

            f(dec_self, *args, **kwargs)
            return _executed

        _tag_function(decorated_with_option)
        _tag_function(decorated_without_option)

        if self.option is None:
            return update_wrapper(decorated_without_option, f)
        else:
            return update_wrapper(decorated_with_option, f)


class Workflow(object):
    """Arbitrary worflow support structure"""

    def __init__(self, short_circuit=True, debug=False, options=None,
                 **kwargs):
        """Build thy self

        short_circuit : if True, enables ignoring function groups when a given
            item has failed
        debug : Enable debug mode
        options : runtime options, {'option':values}
        kwargs : Additional arguments will be added to self

        All workflow methods (i.e., those starting with "wf_") must be
        decorated by either "no_requirements" or "requires". This ensures that
        the methods support the automatic workflow determination mechanism.
        """
        if options is None:
            self.options = {}
        else:
            self.options = options

        ### collections.Counter instead?
        self.stats = defaultdict(int)
        self.short_circuit = short_circuit
        self.failed = False
        self.debug = debug

        if self.debug:
            self.debug_trace = []

        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                raise AttributeError("%s exists in self!" % k)
            setattr(self, k, v)

        for f in self._all_wf_methods():
            if not hasattr(f, '_workflow_tag'):
                raise AttributeError("%s isn't a wf method!" % f.__name__)

        self._allocate_final_state()
        self._setup_debug()

    def _allocate_final_state(self):
        """Setup final_state, must be implemented by subclasses"""
        raise NotImplementedError("Must implement this method")

    def _setup_debug(self):
        """Wrap all methods with debug trace support"""
        if not self.debug:
            return

        # ignore all objects in the baseclass
        ignore = set(dir(Workflow))

        for attrname in dir(self):
            if attrname.startswith('__'):
                continue

            if attrname in ignore:
                continue

            attr = getattr(self, attrname)

            if isinstance(attr, MethodType):
                setattr(self, attrname, _debug_trace_wrapper(self, attr))

    def _all_wf_methods(self, default_priority=0):
        """Get all workflow methods

        Methods are sorted by priority
        """
        methods = [getattr(self, f) for f in dir(self) if f.startswith('wf_')]
        key = lambda x: getattr(x, 'priority', default_priority)
        methods_sorted = sorted(methods, key=key, reverse=True)

        if not self.debug:
            methods_sorted.pop(0)

        return methods_sorted

    def _get_workflow(self, it):
        """Get the methods executed, sorted by priority"""
        # save state
        shortcircuit_state = self.short_circuit
        self.short_circuit = False
        stats = self.stats.copy()

        peek = it.next()
        executed = [f for f in self._all_wf_methods() if f(peek) is _executed]

        # restore state
        self.short_circuit = shortcircuit_state
        self.stats = stats
        generator_reset = chain([peek], it)

        return generator_reset, executed

    @priority(priority.highest)
    @no_requirements
    def wf_setup_debug_trace(self, item):
        self.debug_trace = []

    def __call__(self, it, success_callback=None, fail_callback=None):
        """Operate on all the data

        it : an iterator
        success_callback : method to call on a successful item prior to
            yielding
        fail_callback : method to call on a failed item prior to yielding
        """
        if success_callback is None:
            success_callback = lambda x: x.final_state

        it, workflow = self._get_workflow(it)

        for item in it:
            self.failed = False

            for f in workflow:
                f(item)

            if self.failed:
                if fail_callback is not None:
                    yield fail_callback(self)
            else:
                yield success_callback(self)
