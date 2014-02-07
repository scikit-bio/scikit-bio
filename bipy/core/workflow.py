#!/usr/bin/env python

"""Perform multiple method calls, determined at runtime, on independent items

Construct arbitrarily complex workflows in which the specific methods run are
determined at runtime. These methods are applied to items that are assumed to
be independent.

As an example:

class MyWorkflow(Workflow):
    def _allocate_final_state(self):
        self.FinalState = None

    def _sanity_check(self):
        pass

    @priority(100)
    @no_requirements
    def wf_mul(self, item):
        self.FinalState *= item

    @priority(10)
    @requires(Option='add_value')
    def wf_add(self, item):
        self.FinalState += item

    @requires(Option='sub_value', Values=[1,5,10])
    def wf_sub(self, item):
        self.FinalState -= item
        self.FinalState -= self.Options['sub_value']

    @priority(1000)
    @requires(IsValid=False)
    def wf_init(self, item):
        self.FinalState = item

# (i * i) + i - i - 5
wf = MyWorkflow(Options={'add_value':None, 'sub_value':5})
gen = (i for i in range(10))
for i in wf(gen):
    print i

# (i * i) - i - 10
wf = MyWorkflow(Options={'sub_value':10})
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
        if not hasattr(obj, 'DebugTrace'):
            raise AttributeError("%s doesn't have DebugTrace!" % obj.__class__)

        obj.DebugTrace.append(f.__name__)
        return f(self, *args, **kwargs)

    return update_wrapper(wrapped, f)


def _tag_function(f):
    """Tag, you're it"""
    setattr(f, '_workflow_tag', None)


class priority(object):
    """Decorate a function priority"""
    def __init__(self, Priority):
        self.Priority = Priority

    def __call__(self, f):
        f.Priority = self.Priority
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
    def __init__(self, IsValid=True, Option=None, Values=anything,
                 ValidData=None):
        """
        IsValid : execute the function if self.Failed is False
        Option : a required option
        Values : required values associated with an option
        ValidData : data level requirements, this must be a function with the
            following signature: f(*args, **kwargs) returning True. NOTE: if
            ValidData returns False on the first item evaluated, the decorated
            function will be removed from the remaining workflow.
        """
        # self here is the requires object
        self.IsValid = IsValid
        self.Option = Option
        self.ValidData = ValidData

        if Values is anything:
            self.Values = anything
        elif Values is not_none:
            self.Values = not_none
        elif not isinstance(Values, set):
            if isinstance(Values, str):
                self.Values = Values
            elif isinstance(Values, Iterable):
                self.Values = set(Values)
            else:
                self.Values = set([Values])
        else:
            self.Values = Values

    def doShortCircuit(self, wrapped):
        return self.IsValid and (wrapped.Failed and wrapped.ShortCircuit)

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
            if self.doShortCircuit(dec_self):
                return

            if self.ValidData is not None:
                if not self.ValidData(*args, **kwargs):
                    return

            s_opt = self.Option
            ds_opts = dec_self.Options

            # if the option exists in the Workflow
            if s_opt in ds_opts:
                v = ds_opts[s_opt]

                # if the value just needs to be not None
                if self.Values is not_none:
                    if v is not None:
                        f(dec_self, *args, **kwargs)
                        return _executed

                # otherwise make sure the value is acceptable
                elif ds_opts[s_opt] in self.Values:
                    f(dec_self, *args, **kwargs)
                    return _executed

        def decorated_without_option(dec_self, *args, **kwargs):
            """A decorated function that does not have an option to validate

            dec_self : this is "self" for the decorated function
            """
            if self.doShortCircuit(dec_self):
                return

            if self.ValidData is not None:
                if not self.ValidData(*args, **kwargs):
                    return

            f(dec_self, *args, **kwargs)
            return _executed

        _tag_function(decorated_with_option)
        _tag_function(decorated_without_option)

        if self.Option is None:
            return update_wrapper(decorated_without_option, f)
        else:
            return update_wrapper(decorated_with_option, f)


class Workflow(object):
    """Arbitrary worflow support structure"""

    def __init__(self, ShortCircuit=True, Debug=False, Options=None, **kwargs):
        """Build thy self

        ShortCiruit : if True, enables ignoring function groups when a given
            item has failed
        Debug : Enable debug mode
        Options : runtime options, {'option':values}
        kwargs : Additional arguments will be added to self

        All workflow methods (i.e., those starting with "wk_") must be
        decorated by either "no_requirements" or "requires". This ensures that
        the methods support the automatic workflow determination mechanism.
        """
        if Options is None:
            self.Options = {}
        else:
            self.Options = Options

        ### collections.Counter instead?
        self.Stats = defaultdict(int)
        self.ShortCircuit = ShortCircuit
        self.Failed = False
        self.Debug = Debug

        if self.Debug:
            self.DebugTrace = []

        for k, v in kwargs.iteritems():
            if hasattr(self, k):
                raise AttributeError("%s exists in self!" % k)
            setattr(self, k, v)

        for f in self._all_wf_methods():
            if not hasattr(f, '_workflow_tag'):
                raise AttributeError("%s isn't a wf method!" % f.__name__)

        self._allocate_final_state()
        self._sanity_check()
        self._stage_state()
        self._setup_debug()

    def _allocate_final_state(self):
        """Setup FinalState, must be implemented by subclasses"""
        raise NotImplementedError("Must implement this method")

    def _setup_debug(self):
        """Wrap all methods with debug trace support"""
        if not self.Debug:
            return

        _ignore = set(['_get_workflow', '_all_wf_methods', '_sanity_check',
                       '_stage_state'])

        for attrname in dir(self):
            if attrname.startswith('__'):
                continue

            if attrname in _ignore:
                continue

            attr = getattr(self, attrname)

            if isinstance(attr, MethodType):
                setattr(self, attrname, _debug_trace_wrapper(self, attr))

    def _stage_state(self):
        """Stage any additional data necessary for the workflow

        This does not need to be overloaded
        """
        pass

    def _sanity_check(self):
        """Perform a sanity check on self"""
        raise NotImplementedError("Must implement a sanity check!")

    def _all_wf_methods(self, default_priority=0):
        """Get all workflow methods

        Methods are sorted by priority
        """
        methods = [getattr(self, f) for f in dir(self) if f.startswith('wf_')]
        key = lambda x: getattr(x, 'Priority', default_priority)
        methods_sorted = sorted(methods, key=key, reverse=True)

        if methods_sorted[0] != self.wf_SETUP_DEBUG_TRACE:
            name = methods_sorted[0].__name__
            debug_prio = self.wf_SETUP_DEBUG_TRACE.Priority

            raise AttributeError("Method %s has a higher priority than the "
                                 "debug trace method. Please set its priority "
                                 "below %d." % (name, debug_prio))

        if not self.Debug:
            methods_sorted.pop(0)

        return methods_sorted

    def _get_workflow(self, it):
        """Get the methods executed, sorted by priority"""
        # save state
        shortcircuit_state = self.ShortCircuit
        self.ShortCircuit = False
        stats = self.Stats.copy()

        peek = it.next()
        executed = [f for f in self._all_wf_methods() if f(peek) is _executed]

        # restore state
        self.ShortCircuit = shortcircuit_state
        self.Stats = stats
        generator_reset = chain([peek], it)

        return generator_reset, executed

    @priority(99999999)
    @no_requirements
    def wf_SETUP_DEBUG_TRACE(self, item):
        self.DebugTrace = []

    def __call__(self, it, success_callback=None, fail_callback=None):
        """Operate on all the data

        it : an iterator
        success_callback : method to call on a successful item prior to
            yielding
        fail_callback : method to call on a failed item prior to yielding
        """
        if success_callback is None:
            success_callback = lambda x: x.FinalState

        it, workflow = self._get_workflow(it)

        for item in it:
            self.Failed = False

            for f in workflow:
                f(item)

            if self.Failed:
                if fail_callback is not None:
                    yield fail_callback(self)
            else:
                yield success_callback(self)
