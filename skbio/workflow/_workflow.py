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
