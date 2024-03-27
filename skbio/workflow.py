r"""Workflow construction (:mod:`skbio.workflow`)
=============================================

.. currentmodule:: skbio.workflow

Construct arbitrarily complex workflows in which the specific methods run are
determined at runtime. This module supports short circuiting a workflow if an
item fails, supports ordering methods, callbacks for processed items, and
deciding what methods are executed based on state or runtime options.

Classes
-------

.. autosummary::
    :toctree: generated/

    Workflow


Decorators
----------

.. autosummary::
    :toctree: generated/

    requires
    method


Tutorial
--------

>>> from skbio.workflow import Workflow

As an example of the :class:`Workflow` object, let's construct a sequence processor
that will filter sequences that are < 10 nucleotides, reverse the sequence
if the runtime options indicate to, and truncate if a specific nucleotide
pattern is observed. The :class:`Workflow` object will only short circuit, and
evaluate requirements on methods decorated by ``method``. Developers are free
to define as many methods as they'd like within the object definition, and
which can be called from workflow methods, but they will not be subjected
directly to workflow checks.

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

An instance of a ``Workflow`` must be passed a ``state`` object and any runtime
options. There are a few other useful parameters that can be specfied but are
out of scope for the purposes of this example. We also do not need to provide
a state object as our ``initialize_state`` method overrides ``self.state``.
Now, let's create the instance.

>>> wf = SequenceProcessor(state=None, options={'reverse=': False})

To run items through the ``SequenceProcessor``, we need to pass in an
iterable. So, lets create a ``list`` of sequences.

>>> seqs = ['AAAAAAATTTTTTT', 'ATAGACC', 'AATTGCCGGAC', 'ATATGAACAAA']

Before we run these sequences through, we're going to also define callbacks
that are applied to the result of an single pass through the ``Workflow``.
Callbacks are optional -- by default, a success will simply yield the state
member variable while failures are ignored -- but, depending on your workflow,
it can be useful to handle failures or potentially do something fun and
exciting on success.

>>> def success_f(obj):
...     return "SUCCESS: %s" % obj.state
>>> def fail_f(obj):
...     return "FAIL: %s" % obj.state

Now, lets process some data!

>>> for result in wf(seqs, success_callback=success_f, fail_callback=fail_f):
...     print(result)
SUCCESS: AAAAAAATTTTTTT
FAIL: ATAGACC
SUCCESS: CCGGAC
SUCCESS: ATATGAACAAA

A few things of note just happened. First off, none of the sequences were
reversed as the ``SequenceProcessor`` did not have option "reverse"
set to ``True``. Second, you'll notice that the 3rd sequence was truncated,
which is expected as it matched our nucleotide pattern of interest. Finally,
of the sequences we processed, only a single sequence failed.

To assist in constructing workflows, debug information is available but it
must be turned on at instantiation. Let's do that, and while we're at it, let's
go ahead and enable the reversal method. This time through though, were going
to walk through an item at a time so we can examine the debug information.

>>> wf = SequenceProcessor(state=None, options={'reverse':True}, debug=True)
>>> gen = wf(seqs, fail_callback=lambda x: x.state)
>>> next(gen)
'TTTTTTTAAAAAAA'
>>> wf.failed
False
>>> sorted(wf.debug_trace)
[('check_length', 0), ('reverse', 2)]

The ``debug_trace`` specifies the methods executed, and the order of their
execution where closer to zero indicates earlier in the execution order. Gaps
indicate there was a method evaluated but not executed. Each of the items in
the ``debug_trace`` is a key into a few other ``dict`` of debug information
which we'll discuss in a moment. Did you see that the sequence was reversed
this time through the workflow?

Now, let's take a look at the next item, which on our prior run through the
workflow was a failed item.

>>> next(gen)
'ATAGACC'
>>> wf.failed
True
>>> sorted(wf.debug_trace)
[('check_length', 0)]

What we can see is that the failed sequence only executed the check_length
method. Since the sequence didn't pass our length filter of 10 nucleotides,
it was marked as failed within the ``check_length`` method. As a result, none
of the other methods were evaluated (note: this short circuiting behavior can
be disabled if desired).

This third item previously matched our nucleotide pattern of interest for
truncation. Let's see what that looks like in the debug output.

>>> next(gen)
'CAGGCC'
>>> wf.failed
False
>>> sorted(wf.debug_trace)
[('check_length', 0), ('reverse', 2), ('truncate', 1)]

In this last example, we can see that the ``truncate`` method was executed
prior to the ``reverse`` method and following the ``check_length`` method. This
is as anticipated given the priorities we specified for these methods. Since
the ``truncate`` method is doing something interesting, let's take a closer
look at how the ``state`` is changing. First, we're going to dump out the
state of the workflow prior to the call to ``truncate`` and then we're going
to dump out the ``state`` following the call to ``truncate``, which will allow
us to rapidly what is going on.

>>> wf.debug_pre_state[('truncate', 1)]
'AATTGCCGGAC'
>>> wf.debug_post_state[('truncate', 1)]
'CCGGAC'

As we expect, we have our original sequence going into ``truncate``, and
following the application of ``truncate``, our sequence is missing our
nucleotide pattern of interest. Awesome, right?

There is one final piece of debug output, ``wf.debug_runtime``, which can
be useful when diagnosing the amount of time required for individual methods
on a particular piece of state (as opposed to the aggregate as provided by
cProfile).

Three final components of the workflow that are quite handy are objects that
allow you to indicate ``anything`` as an option value, anything that is
``not_none``, and a mechanism to define a range of valid values.

>>> from skbio.workflow import not_none, anything
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


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from copy import deepcopy
from time import time
from functools import update_wrapper
from collections.abc import Iterable
from types import MethodType


class NotExecuted:
    """Helper object to track if a method was executed."""

    def __init__(self):
        """Construct all the necessary attributes for the NotExecuted object."""
        self.msg = None

    def __call__(self, msg):
        """Update message and return self."""
        self.msg = msg
        return self


_not_executed = NotExecuted()


class Exists:
    """Stub object to assist with ``requires`` when a value exists."""

    def __contains__(self, item):
        """Check if a value exists."""
        return True


anything = Exists()  # external, for when a value can be anything


class NotNone:
    """Check for non-None values."""

    def __contains__(self, item):
        """Check if item is not None."""
        if item is None:
            return False
        else:
            return True


not_none = NotNone()


class Workflow:
    """Arbitrary workflow support structure.

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

    """

    def __init__(self, state, short_circuit=True, debug=False, options=None, **kwargs):
        r"""Build thy workflow of self."""
        if options is None:
            self.options = {}
        else:
            self.options = options

        self.short_circuit = short_circuit
        self.failed = False
        self.debug = debug
        self.state = state
        self.iter_ = None

        for k, v in kwargs.items():
            if hasattr(self, k):
                raise AttributeError("'%s' already exists in self." % k)
            setattr(self, k, v)

        if self.debug:
            self._setup_debug()

    def initialize_state(self, item):
        """Initialize state.

        This method is called first prior to any other defined workflow method
        with the exception of _setup_debug_trace if self.debug is True

        Parameters
        ----------
        item : anything
            Workflow dependent

        """
        raise NotImplementedError("Must implement this method")

    def _setup_debug(self):
        """Wrap all methods with debug trace support."""
        # ignore all members of the baseclass
        ignore = set(dir(Workflow))

        for attrname in dir(self):
            if attrname in ignore:
                continue

            attr = getattr(self, attrname)

            if isinstance(attr, MethodType):
                setattr(self, attrname, self._debug_trace_wrapper(attr))

    def _all_wf_methods(self):
        """Get all workflow methods.

        Methods are sorted by priority
        """
        methods = []
        for item in dir(self):
            obj = getattr(self, item)
            if hasattr(obj, "priority"):
                methods.append(obj)

        def key(x):
            return getattr(x, "priority")

        methods_sorted = sorted(methods, key=key, reverse=True)

        if self.debug:
            methods_sorted.insert(0, self._setup_debug_trace)

        return methods_sorted

    def _setup_debug_trace(self):
        """Set up a trace.

        The trace is per item iterated over by the workflow. Information about
        each method executed is tracked and keyed by::

            (function name, order of execution)

        Order of execution starts from zero. Multiple calls to the same
        function are independent in the trace.

        The following information is tracked::

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
        """Operate on all the data.

        This is the processing engine of the workflow. Callbacks are executed
        following applying all workflow methods to an item from ``iter_``
        (unless ``short_cicruit=True`` in which case method execution for an
        item is stopped if ``failed=True``). Callbacks are provided ``self``
        which allows them to examine any aspect of the workflow.

        Parameters
        ----------
        iter_ : iterator
            The iterator containing the data to be processed.
        success_callback : method, optional
            Method to call on a successful item prior to yielding. By default,
            ``self.state`` is yielded.
        fail_callback : method, optional
            Method to call on a failed item prior to yielding. By default, failures
            are ignored.

        """
        if success_callback is None:

            def success_callback(x):
                return x.state

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
        """Trace a function call."""

        def wrapped():
            """Track debug information about a method execution."""
            if not hasattr(self, "debug_trace"):
                raise AttributeError("%s doesn't have debug_trace." % self.__class__)

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


class method:
    """Decorate a function to indicate it is a workflow method.

    Parameters
    ----------
    priority : int
        Specify a priority for the method, the higher the value the higher
        the priority. Priorities are relative to a given workflow

    """

    highest_priority = sys.maxsize

    def __init__(self, priority=0):
        """Construct all the necessary attributes for the method object."""
        self.priority = priority

    def __call__(self, func):
        """Decorate function with specified priority."""
        func.priority = self.priority
        return func


class requires:
    """Decorator that executes a function if requirements are met.

    Parameters
    ----------
    option : any Hashable object
        An option that is required for the decorated method to execute.
        This option will be looked up within the containing ``Workflow``s'
        ``options``.
    values : object
        A required value. This defaults to ``anything`` indicating that
        the only requirement is that the ``option`` exists. It can be
        useful to specify ``not_none`` which indicates that the
        requirement is satisfied if the ``option`` exists and it holds
        a value that is not ``None``. Values also supports iterables
        or singular values.
    state : Function
        A requirement on workflow state. This must be a function that
        accepts a single argument, and returns ``True`` to indicate
        the requirement is satisfied, or ``False`` to indicate the
        requirement is not satisfied. This method will be passed the
        containing ``Workflow``s' ``state`` member variable.

    """

    def __init__(self, option=None, values=anything, state=None):
        """Construct all the necessary attributes for the requires object."""
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
        """Wrap a function.

        func : the function to wrap
        """

        def decorated(dec_self):
            """Execute a decorated function that has requirements.

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
