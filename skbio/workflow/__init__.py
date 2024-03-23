r"""Constructing workflows (:mod:`skbio.workflow`)
==============================================

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

Examples
--------
>>> from skbio.workflow import Workflow

As an example of the ``Workflow`` object, let's construct a sequence processor
that will filter sequences that are < 10 nucleotides, reverse the sequence
if the runtime options indicate to, and truncate if a specific nucleotide
pattern is observed. The ``Workflow`` object will only short circuit, and
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

from ._workflow import Workflow, requires, method, not_none, anything


__all__ = ["Workflow", "requires", "method", "not_none", "anything"]
