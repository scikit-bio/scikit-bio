Supporting Python 2 and Python 3
################################

skbio simultaneously supports Python 2.7 and 3.3+ by writing code that
works unchanged in both major versions.

As a compatibility layer, we're using the `future <http://python-future.org/>`_
and `six <https://pypi.python.org/pypi/six>`_ projects. future "allows you to
use a single, clean Python 3.x-compatible codebase to support both Python 2 and
Python 3 with minimal overhead". It includes functionality from "six, IPython,
Jinja2, Django, and Pandas". Recent versions of the future project stopped
bundling the six library, so we also directly depend on six (e.g., for StringIO
compatibility).

So far, these notes are based on issues that have appeared when porting
skbio, so it is not a complete guide. Refer to the `official porting
guide <https://docs.python.org/3/howto/pyporting.html>`_ and the
`python-future docs <http://python-future.org/>`_ for more
information.

Importing __future__
====================

For consistency across versions, every Python file should start with
the following imports::

  # ----------------------------------------------------------------------------
  # Copyright (c) 2013--, scikit-bio development team.
  #
  # Distributed under the terms of the Modified BSD License.
  #
  # The full license is in the file COPYING.txt, distributed with this software.
  # ----------------------------------------------------------------------------

  from __future__ import absolute_import, division, print_function

Iterators
=========

Builtins
--------

Builtin iterators in Python 2 usually return lists, and have an
alternative that returns an iterator (i.e., `range` and `xrange`,
`items` and `iteritems`). In Python 3, only the iterator version
exists but it uses the list-returning name (i.e., `range` and
`items`).

When iterating over the resulting object, the recommended approach
depends on efficiency concerns:

- If iteration only happens over a few items, you can use the
  functions that exist both in Python 2 and Python 3.

- If the number of iterations can be large and efficiency is
  important, use the future package.

+--------------------+----------------------------+--------------------+
|Small # of          |Efficient versions          |Notes               |
|iterations (returns |(always iterators)          |                    |
|lists in py2,       |                            |                    |
|iterators in py3)   |                            |                    |
+--------------------+----------------------------+--------------------+
|`zip`               |`future.builtins.zip`       |                    |
+--------------------+----------------------------+--------------------+
|`range`             |`future.builtins.range`     |                    |
+--------------------+----------------------------+--------------------+
|`map`               |`future.builtins.map`       |Prefer lists        |
|                    |                            |comprehensions or   |
|                    |                            |for loops in        |
|                    |                            |general. Avoid      |
|                    |                            |calling functions   |
|                    |                            |that cause side     |
|                    |                            |effects when using  |
|                    |                            |map. Gotcha: Py3's  |
|                    |                            |`map` stops when the|
|                    |                            |shortest iterable is|
|                    |                            |exhausted, but Py2's|
|                    |                            |pads them with      |
|                    |                            |`None` till the     |
|                    |                            |longest iterable is |
|                    |                            |exhausted.          |
|                    |                            |                    |
+--------------------+----------------------------+--------------------+
|`filter`            |`future.builtins.filter`    |                    |
|                    |                            |                    |
|                    |                            |                    |
+--------------------+----------------------------+--------------------+
|`functools.reduce`  |`functools.reduce`          |Avoid using the     |
|                    |                            |global reduce       |
|                    |                            |available in Py2 (it|
|                    |                            |is the same as the  |
|                    |                            |`functools` one)    |
+--------------------+----------------------------+--------------------+
|`d.items()`         |`future.utils.viewitems(d)` |Efficient iteration |
|                    |                            |over d *and*        |
|                    |                            |set-like behaviour  |
+--------------------+----------------------------+--------------------+
|`d.values()`        |`future.utils.viewvalues(d)`|Efficient iteration |
|                    |                            |over d *and*        |
|                    |                            |set-like behaviour  |
+--------------------+----------------------------+--------------------+
|`d.keys()`          |`future.utils.viewkeys(d)`  |Hardly ever needed, |
|                    |                            |as iterating over a |
|                    |                            |dictionary yields   |
|                    |                            |keys (thus sorted(d)|
|                    |                            |returns the sorted  |
|                    |                            |keys).              |
+--------------------+----------------------------+--------------------+


When not directly iterating over an iterator, don't write code that
relies on list-like behaviour: you may need to cast it explicitly. The
following snippets show some possible issues::

    a = zip(...)
    b = zip(...)
    c = a + b  # succeeds in Py2 (list concatenation), TypeError in Py3

::

    s = map(int, range(2))
    1 in s  # True (membership testing in a list is an O(n) bad idea)
    0 in s  # True in Py2, False in Py3

In Py2, `s` is a list, so clearly `(1 in [0, 1]) == True` and `(0 in
[0, 1]) == True`. In Py3, `s` is an iterator and the items it yields
are discarded. Let's see an example with a generator to try and make
it more clear::

    >>> s = ((i, print(i)) for i in [0, 1, 2])  # print will let us see the iteration
    >>> (1, None) in s  # Starts iterating over s...
    0
    1                   # ...till it finds (1, None)
    True
    >>> (0, None) in s  # Continues iterating over s
    2                   # s is exhausted
    False               # but (0, None) isn't there


Advancing an iterator
---------------------

Always use the next function, which is available from Python 2.6
onwards. Never call the next method, which doesn't exist in Py3.

Implementing new iterators
--------------------------

Implement the `__next__` special method, like in Py3, and decorate the
class::

    from future.utils import implements_iterator

    @implements_iterator
    class ParameterIterBase(object):
    def __next__(self):
        return next(self._generator)

It is also possible to subclass from `future.builtins.object`. In this
case, no decorator is needed.

Changes in the standard library
===============================

To deal with modules that live under a different place, future
provides a context manager::

    # Example from future's documentation
    from future import standard_library

    with standard_library.hooks():
        from http.client import HttpConnection
        from itertools import filterfalse
        import html.parser
        import queue

StringIO and BytesIO
--------------------

In Py2 there are three flavours of StringIO: a pure Python module
(StringIO), an accelerated version (cStringIO), and another one in the
io module. They all behave in a slightly different way, with different
memory and performance characteristics. So far, we're using::

    from six import StringIO

It refers to `io.StringIO` in Py3, and `StringIO.StringIO` in Py2.

If you need a binary file-like object (see the Text vs bytes section),
use `six.BytesIO`, which refers to `io.BytesIO` in Py3, and `StringIO.StringIO`
in Py2.

Text vs bytes
=============

This is a fundamental change between Py2 and Py3. It is very important
to always distinguish text from bytes.

String literals that are to be treated as bytes need the `b`
prefix. String literals that are text need either the `u` prefix or
`from __future__ import unicode_literals` at the top.

A brief introduction: Unicode, UTF-8, ASCII...
----------------------------------------------

A string can be seen as a sequence of characters. According to the
Unicode standard, each character is represented by a code point (a
number). For example, character `ñ` is represented by the Unicode code
point `U+00F1`. Code points are still abstract and can be stored in a
number of ways, including even little or big endian formats. There are
many encodings that map code points to byte values (encode) and back
(decode). Three important ones are ASCII, UTF-8 and latin-1:

- ASCII is a 7 bit encoding that can handle a very limited range of
  Unicode code points (not even the one corresponding to character
  `ñ`).

- UTF-8 is an encoding that can represent every Unicode character. It
  is ASCII-compatible because code points that can also be represented
  by ASCII are mapped to the same byte value by UTF-8 and ASCII. `ñ`
  is represented by the byte sequence `\xC3\xB1`.

- latin-1 is an ASCII-compatible 8 bit encoding that maps the first
  256 Unicode code points to their byte values. That is, the Unicode
  code point `U+00F1` (character `ñ`) is directly encoded as `0xF1` in
  latin-1. The Py2 `str` type loosely worked by assuming everything
  was encoded in latin-1.


Text processing
---------------

    There Ain't No Such Thing As Plain Text.  -- Joel Spolsky, `The
    Absolute Minimum Every Software Developer Absolutely, Positively
    Must Know About Unicode and Character Sets (No Excuses!)
    <http://www.joelonsoftware.com/articles/Unicode.html>`_, 2003.

After going through Nick Coghlan's `"Processing Text Files in Python
3"
<https://ncoghlan_devs-python-notes.readthedocs.org/en/latest/python3/text_file_processing.html>`_
I think the way forward is to process ASCII-like files (fasta, fastq)
as binary files, and decode to strings some parts, if necessary. This
is faster than processing them as text files, especially in Py3. In
fact, it seems (from functions like `_phred_to_ascii*`) that these
formats are in fact mixed binary and ASCII, which I think puts us in
the same place as people dealing with `network protocols
<https://ncoghlan_devs-python-notes.readthedocs.org/en/latest/python3/binary_protocols.html>`_:
it's more cumbersome to do in Py3, especially before Python 3.5
arrives, which will `reintroduce binary string interpolation
<http://legacy.python.org/dev/peps/pep-0460/>`_).

Gotchas
-------

Comparing bytes and text strings always returns `False` in Python 3
(as they're incompatible types, and comparisons are required to
succeed by the language)::

    >>> b'GATCAT' == 'GATCAT'
    False

Calling `str` on a bytes instance returns a string with the `b` prefix
and quotes, which will give unexpected results when using string
formatting::

    >>> "Sequence {}".format(b'GATCAT')
    "Sequence b'GATCAT'"

If you actually want to construct a text string, bytes objects need to
be *decoded* into text. For example::

    >>> "Sequence {}".format(b'GATCAT'.decode('utf-8'))

If you want to efficiently construct a byte string, the most
convenient way may be to call `b''.join(iterable of byte strings)`,
though there are other options like using `io.BytesIO` or
`bytearray`. For a very small number of byte strings, it may be OK to
use the `+` operator.

Run python with the `-b` flag to detect these two bug-prone usages,
and `-bb` to turn them into exceptions.

Instance checking: basestring, str, unicode, bytes, long, int
=============================================================

Strings
-------

When testing if a variable is a string use
`six.string_types`. It refers to `basestring` in Py2 and `str` in Py3.
`binary_type` and `text_type` are also available.

Numbers
-------

The `long` type no longer exists in Py3. To test if a number is an
integer (`int` or `long` in Py2, `int` in Py3), compare it to
the abstract base class `Integral`::

    from numbers import Integral
    isinstance(quality, Integral)

Implementing comparisons
========================

If the class you're defining has a `total ordering
<http://en.wikipedia.org/wiki/Total_order>`_, either use
`functools.total_ordering
<https://docs.python.org/2.7/library/functools.html#functools.total_ordering>`_
or implement all rich comparison methods if comparison performance is
a bottleneck. Don't implement `__cmp__`, which was removed in Py3.

However, usually only equality is important and you should only define
`__eq__`. While compatibility with Py2 is kept, `__ne__` needs to be
implemented too::

    def __ne__(self, other):
        """Required in Py2."""
        return not self == other

Otherwise, using the operator `!=` will lead to unexpected results in
Py2 because it will compare identity, not equality::

    class Foo(object):
        def __eq__(self, other):
            return True

    print(Foo() != Foo())

That prints `True` in Py2 (because each instance has a different `id`)
but prints `False` in Py3 (the opposite of what `__eq__` returns,
which is the desired behaviour).

Always test that both `==` and `!=` are behaving correctly, e.g.::

    def test_eq(self):
        gc_1 = GeneticCode(self.sgc)
        gc_2 = GeneticCode(self.sgc)
        self.assertEqual(gc_1, gc_2)

    def test_ne(self):
        gc_1 = GeneticCode(self.sgc)
        gc_2 = GeneticCode(self.sgc)
        # Explicitly using !=
        self.assertFalse(gc_1 != gc_2)

Other modules
=============

Numpy
-----

Try to avoid setting dtypes to a string (i.e., use `dtype=np.float64`
instead of `dtype='float'`, etc). It is may be safe, but some warnings
were raised when running Python with the `-b` flag. Also, field names
in structured dtypes need to be bytes (`str` type) in Py2, but text
(`str` type) in Py3 (`issue #2407
<https://github.com/numpy/numpy/issues/2407>`_).

Testing
=======

`unittest.assertEquals` is deprecated. Use `unittest.assertEqual`
instead. The complete list of deprecated testing methods is `here
<https://docs.python.org/3.4/library/unittest.html#deprecated-aliases>`_
