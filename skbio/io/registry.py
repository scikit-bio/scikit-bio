r"""I/O Registry (:mod:`skbio.io.registry`)
=======================================

.. currentmodule:: skbio.io.registry

Classes
-------

.. autosummary::
   :toctree:

   IORegistry
   Format

Functions
---------

.. autosummary::
   :toctree:

   create_format

Exceptions
----------

.. autosummary::

   DuplicateRegistrationError
   InvalidRegistrationError


Creating a new format for scikit-bio
------------------------------------
scikit-bio makes it simple to add new file formats to its I/O registry.
scikit-bio maintains a singleton of the :class:`IORegistry` class called
`io_registry`. This is where all scikit-bio file formats are registered. One
could also instantiate their own :class:`IORegistry`, but that is not the focus
of this tutorial.

The first step to creating a new format is to add a submodule in
`skbio/io/format/` named after the file format you are implementing.
For example, if the format you are implementing is called `myformat` then you
would create a file called `skbio/io/format/myformat.py`.

The next step is to import the :func:`create_format` factory from
:mod:`skbio.io`. This will allow you to create a new :class:`Format` object
that `io_registry` will know about.

Ideally you should name the result of :func:`create_format` as your file name.
For example:

.. code-block:: python

   from skbio.io import create_format

   myformat = create_format('myformat')

The `myformat` object is what we will use to register our new functionality.
At this point you should evaulate whether your format is binary or text.
If your format is binary, your :func:`create_format` call should look like
this:

.. code-block:: python

   myformat = create_format('myformat', encoding='binary')

Alternatively if your format is text and has a specific encoding or newline
handling you can also specify that:

.. code-block:: python

   myformat = create_format('myformat', encoding='ascii', newline='\n')

This will ensure that our registry will open files with a default encoding of
`'ascii'` for `'myformat'` and expect all newlines to be `'\n'` characters.

Having worked out these details, we are ready to register the actual
functionality of our format (e.g., sniffer, readers, and writers).

To create a sniffer simply decorate the following onto your sniffer function:

.. code-block:: python

   @myformat.sniffer()
   def _myformat_sniffer(fh):
       # do something with `fh` to determine the membership of the file

For futher details on sniffer functions see :func:`Format.sniffer`.

Creating a reader is very similar, but has one difference:

.. code-block:: python

   @myformat.reader(SomeSkbioClass)
   def _myformat_to_some_skbio_class(fh, kwarg1='default', extra=FileSentinel):
       # parse `fh` and return a SomeSkbioClass instance here
       # `extra` will also be an open filehandle if provided else None

Here we bound a function to a specific class. We also demonstrated using
our FileSentinel object to indicate to the registry that this reader can take
auxilary files that should be handled in the same way as the primary file.
For futher details on reader functions see :func:`Format.reader`.

Creating a writer is about the same:

.. code-block:: python

   @myformat.writer(SomeSkbioClass)
   def _some_skbio_class_to_myformat(obj, fh, kwarg1='whatever',
                                     extra=FileSentinel):
       # write the contents of `obj` into `fh` and whatever else into `extra`
       # do not return anything, it will be ignored

This is exactly the same as the `reader` above just in reverse, we also
receive the object we are writing as the first parameter instead of the file
(which is the second one). For further details on writer functions see
:func:`Format.writer`.

.. note:: When raising errors in readers and writers, the error should be a
   subclass of ``FileFormatError`` specific to your new format.

Once you are satisfied with the functionality, you will need to ensure that
`skbio/io/__init__.py` contains an import of your new submodule so the
decorators are executed. Add the function
``import_module('skbio.io.format.myformat')`` with your module name to the
existing list.

.. note:: Because scikit-bio handles all of the I/O boilerplate, you only need
   to unit-test the actual business logic of your `readers`, `writers`, and
   `sniffers`.

Reserved Keyword Arguments
--------------------------
The following keyword args may not be used when defining new `readers` or
`writers` as they already have special meaning to the registry system:

- `format`
- `into`
- `verify`
- `mode`
- `encoding`
- `errors`
- `newline`
- `compression`
- `compresslevel`

The following are not yet used but should be avoided as well:

- `auth`
- `user`
- `password`
- `buffering`
- `buffer_size`
- `closefd`
- `exclusive`
- `append`

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
import types
import traceback
import itertools
import inspect
from functools import wraps

from ._exception import DuplicateRegistrationError, InvalidRegistrationError
from . import (
    UnrecognizedFormatError,
    ArgumentOverrideWarning,
    FormatIdentificationWarning,
)
from .util import _resolve_file, open_file, open_files, _d as _open_kwargs
from skbio.util._misc import make_sentinel, find_sentinels
from skbio.util._decorator import classonlymethod

FileSentinel = make_sentinel("FileSentinel")


class IORegistry:
    """Create a registry of formats and implementations which map to classes."""

    def __init__(self):
        """Initialize registry mapping formats and implementations to classes.

        This seperation of binary and text formats is useful because there
        are many situations where we may have recieved a text-file. When this
        happens, the binary data fundamentally does not exist. We could
        assume encoding should be interpreted in reverse, however this misses
        the bigger point: why would the user ever want text to be treated as
        binary? They already went through the effort to hand us text.
        Therefore, during format resolution, we should skip the binary
        formats if they are irrelevant. (They are incompatible with such a
        filehandle anyways.)
        """
        self._binary_formats = {}
        self._text_formats = {}
        self._lookups = (self._binary_formats, self._text_formats)

    def create_format(self, *args, **kwargs):
        """Create new file formats.

        A simple factory for creating new file formats.
        This will automatically register the format with this regsistry.
        All arguments are passed through to the Format constructor.

        Returns
        -------
        Format
            A new format that is registered with the registry.

        """
        format = Format(*args, **kwargs)
        self.add_format(format)
        return format

    def add_format(self, format_object):
        """Add a format to the registry.

        Parameters
        ----------
        format_object : Format
            The format to add to the registry.

        """
        # See comment in the constructor for an explanation for why this split
        # occurs.
        name = format_object.name
        if name in self._binary_formats or name in self._text_formats:
            raise DuplicateRegistrationError(
                "A format already exists with" " that name: %s" % name
            )

        if format_object.is_binary_format:
            self._binary_formats[name] = format_object
        else:
            self._text_formats[name] = format_object

    def get_sniffer(self, format_name):
        """Locate the sniffer for a format.

        Parameters
        ----------
        format_name : str
            The name of the format to lookup.

        Returns
        -------
        function or None
            The sniffer associated with `format_name`

        """
        for lookup in self._lookups:
            if format_name in lookup:
                return lookup[format_name].sniffer_function
        return None

    def get_reader(self, format_name, cls):
        """Locate the reader for a format and class.

        Parameters
        ----------
        format_name : str
            The name of the format to lookup.
        cls : type or None
            The class which the reader will return an instance of. If `cls` is
            None, the reader will return a generator.
            Default is None.

        Returns
        -------
        function or None
            The reader associated with `format_name` and `cls`

        """
        return self._get_rw(format_name, cls, "readers")

    def get_writer(self, format_name, cls):
        """Locate the writer for a format and class.

        Parameters
        ----------
        format_name : str
            The name of the format to lookup.
        cls : type or None
            The class which the writer will expect an instance of. If `cls` is
            None, the writer will expect a generator.
            Default is None.

        Returns
        -------
        function or None
            The writer associated with `format_name` and `cls`

        """
        return self._get_rw(format_name, cls, "writers")

    def _get_rw(self, format_name, cls, lookup_name):
        for lookup in self._lookups:
            if format_name in lookup:
                format_lookup = getattr(lookup[format_name], lookup_name)
                if cls in format_lookup:
                    return format_lookup[cls]
        return None

    def list_read_formats(self, cls):
        """Return a list of available read formats for a given `cls` type.

        Parameters
        ----------
        cls : type
            The class which will be used to determine what read formats exist
            for an instance of `cls`.

        Returns
        -------
        list
            A list of available read formats for an instance of `cls`. List may
            be empty.

        """
        return list(self._iter_rw_formats(cls, "readers"))

    def list_write_formats(self, cls):
        """Return a list of available write formats for a given `cls` type.

        Parameters
        ----------
        cls : type
            The class which will be used to determine what write formats exist
            for an instance of `cls`.

        Returns
        -------
        list
            A list of available write formats for an instance of `cls`. List
            may be empty.

        """
        return list(self._iter_rw_formats(cls, "writers"))

    def _iter_rw_formats(self, cls, lookup_name):
        for lookup in self._lookups:
            for format in lookup.values():
                if cls in getattr(format, lookup_name):
                    yield format.name

    def sniff(self, file, **kwargs):
        r"""Detect the format of a given file and suggest kwargs for reading.

        Parameters
        ----------
        file : openable (filepath, URL, filehandle, etc.)
            The file to sniff. Something that is understood by :func:`skbio.io.open`.
        kwargs : dict, optional
            Keyword arguments will be passed to :func:`skbio.io.open`. ``newline``
            cannot be provided.

        Returns
        -------
        (str, dict)
            The name of the format of the file and any suggested kwargs for
            use with the corresponding reader.

        Raises
        ------
        UnrecognizedFormatError
            This occurs when the format is not 'claimed' by any registered
            sniffer or when the format is ambiguous and has been 'claimed' by
            more than one sniffer.
        TypeError
            If ``newline`` is provided in ``kwargs``.

        """
        if "newline" in kwargs:
            raise TypeError("Cannot provide `newline` keyword argument when sniffing.")

        # By resolving the input here, we have the oppurtunity to reuse the
        # file (which is potentially ephemeral). Each sniffer will also resolve
        # the file, but that call will short-circuit and won't claim
        # responsibility for closing the file. This means that the file
        # should only close after leaving this context. This is also the reason
        # that we have to use SaneTextIOWrapper because each sniffer will
        # wrap the file to produce an appropriate default encoding for their
        # format (if unspecified). This results in the SaneTextIOWrapper being
        # garbage collected (using io.TextIOBase results in close being called
        # on our buffer by the deconstructor which we wanted to share with the
        # next sniffer)
        with _resolve_file(file, mode="r", **kwargs) as (fh, _, is_binary_file):
            # tell may fail noisily if the user provided a TextIOBase or
            # BufferedReader which has already been iterated over (via next()).
            matches = []
            backup = fh.tell()
            if is_binary_file and kwargs.get("encoding", "binary") == "binary":
                matches = self._find_matches(fh, self._binary_formats, **kwargs)

            if kwargs.get("encoding", None) != "binary":
                # We can always turn a binary file into a text file, but the
                # reverse doesn't make sense.
                matches += self._find_matches(fh, self._text_formats, **kwargs)
                fh.seek(backup)
            elif not is_binary_file:
                raise ValueError("Cannot decode text source (%r) as binary." % file)
            # else we are a binary_file and our encoding did not exclude binary
            # so we have already handled that condition

        if len(matches) > 1:
            raise UnrecognizedFormatError(
                "File format for %r is ambiguous,"
                " may be one of: %r" % (file, [m for m, s in matches])
            )
        elif len(matches) == 0:
            raise UnrecognizedFormatError("Could not detect the format of %r" % file)

        return matches[0]

    def _find_matches(self, file, lookup, **kwargs):
        matches = []
        for format in lookup.values():
            if format.sniffer_function is not None:
                is_format, skwargs = format.sniffer_function(file, **kwargs)
                file.seek(0)
                if is_format:
                    matches.append((format.name, skwargs))
        return matches

    def read(self, file, format=None, into=None, verify=True, **kwargs):
        r"""Read a file as certain format into an object.

        Parameters
        ----------
        file : openable (filepath, URL, filehandle, etc.)
            The file to read. Something that is understood by :func:`skbio.io.open`.
        format : str, optional
            The format of the file if known. If None, the format will be inferred from
            the file.
        into : type or None, optional
            The object which will be returned. If None, a generator will be returned.
        verify : bool, optional
            When True, will double check the ``format`` if provided.
        kwargs : dict, optional
            Keyword arguments will be passed to their respective handlers
            (:func:`skbio.io.open` and the reader for ``format``). ``newline`` cannot
            be provided.

        Returns
        -------
        object or generator
            An instance of ``into`` if ``into`` is not None else generator

        Raises
        ------
        ValueError
            Raised when ``format`` and ``into`` are both None.
        TypeError
            If ``newline`` is provided in ``kwargs``.
        UnrecognizedFormatError
            Raised when a reader could not be found for a given ``format`` or the
            format could not be guessed.
        FormatIdentificationWarning
            Raised when ``verify`` is True and the sniffer of a ``format`` did
            not agree that ``file`` is a member of ``format``
        ArgumentOverrideWarning
            Raised when ``verify`` is True and a user-supplied argument is
            overriding the suggestion provided by the sniffer of ``format``.

        """
        if "newline" in kwargs:
            raise TypeError("Cannot provide `newline` keyword argument when reading.")

        # Context managers do not compose well with generators. We have to
        # duplicate the logic so that the file will stay open while yielding.
        # Otherwise the context exits as soon as the generator is returned
        # (making any iteration fail as the file is closed from its
        # perspective).
        if into is None:
            if format is None:
                raise ValueError("`into` and `format` cannot both be None")
            gen = self._read_gen(file, format, into, verify, kwargs)
            # This is done so that any errors occur immediately instead of
            # on the first call from __iter__
            # eta-reduction is possible, but we want to the type to be
            # GeneratorType
            try:
                return (x for x in itertools.chain([next(gen)], gen))
            except StopIteration:
                # If the error was a StopIteration, then we want to return an
                # empty generator as `next(gen)` failed.
                # See #1313 for more info.
                return (x for x in [])
        else:
            return self._read_ret(file, format, into, verify, kwargs)

    def _read_ret(self, file, fmt, into, verify, kwargs):
        io_kwargs = self._find_io_kwargs(kwargs)
        with _resolve_file(file, **io_kwargs) as (file, _, _):
            reader, kwargs = self._init_reader(
                file, fmt, into, verify, kwargs, io_kwargs
            )
            return reader(file, **kwargs)

    def _read_gen(self, file, fmt, into, verify, kwargs):
        io_kwargs = self._find_io_kwargs(kwargs)
        # We needed to get the io_kwargs from kwargs for things like
        # _resolve_file and for verifying a format.
        # kwargs should still retain the contents of io_kwargs because the
        # actual reader will also need them.
        with _resolve_file(file, **io_kwargs) as (file, _, _):
            reader, kwargs = self._init_reader(
                file, fmt, into, verify, kwargs, io_kwargs
            )
            yield from reader(file, **kwargs)

    def _find_io_kwargs(self, kwargs):
        return {k: kwargs[k] for k in _open_kwargs if k in kwargs}

    def _init_reader(self, file, fmt, into, verify, kwargs, io_kwargs):
        skwargs = {}
        if fmt is None:
            fmt, skwargs = self.sniff(file, **io_kwargs)
        elif verify:
            sniffer = self.get_sniffer(fmt)
            if sniffer is not None:
                backup = file.tell()
                is_format, skwargs = sniffer(file, **io_kwargs)
                file.seek(backup)
                if not is_format:
                    warn(
                        "%r does not look like a %s file" % (file, fmt),
                        FormatIdentificationWarning,
                    )

        for key in skwargs:
            if key not in kwargs:
                kwargs[key] = skwargs[key]
            elif kwargs[key] != skwargs[key]:
                warn(
                    "Best guess was: %s=%r, continuing with user"
                    " supplied: %r" % (key, skwargs[key], kwargs[key]),
                    ArgumentOverrideWarning,
                )

        reader = self.get_reader(fmt, into)
        if reader is None:
            possible_intos = [r.__name__ for r in self._get_possible_readers(fmt)]
            message = ""
            if possible_intos:
                message = "Possible values for `into` include: %s" % ", ".join(
                    possible_intos
                )
            into_message = "`into` also not provided." if not into else ""
            raise UnrecognizedFormatError(
                "Cannot read %r from %r, no %s reader found. %s %s"
                % (
                    fmt,
                    file,
                    into.__name__ if into else "generator",
                    into_message,
                    message,
                )
            )
        return reader, kwargs

    def _get_possible_readers(self, fmt):
        for lookup in self._lookups:
            if fmt in lookup:
                return list(lookup[fmt].readers)
        return []

    def write(self, obj, format, into, **kwargs):
        r"""Write an object as certain format into a file.

        Parameters
        ----------
        obj : object
            The object to write as ``format``.
        format : str
            The format to write ``obj`` as.
        into : openable (filepath, URL, filehandle, etc.)
            What to write ``obj`` to. Something that is understood by
            :func:`skbio.io.open`.
        kwargs : dict, optional
            Keyword arguments will be passed to their respective handlers
            (:func:`skbio.io.open` and the writer for ``format``).

        Returns
        -------
        openable (filepath, URL, filehandle, etc.)
            Will pass back the user argument for ``into`` as a convenience.

        Raises
        ------
        UnrecognizedFormatError
            Raised when a writer for writing ``obj`` as ``format`` could not be found.

        """
        # The simplest functionality here.
        cls = None
        if not isinstance(obj, types.GeneratorType):
            cls = obj.__class__
        writer = self.get_writer(format, cls)
        if writer is None:
            raise UnrecognizedFormatError(
                "Cannot write %r into %r, no %s writer found."
                % (format, into, obj.__class__.__name__)
            )

        writer(obj, into, **kwargs)
        return into

    def monkey_patch(self):
        r"""Monkey-patch ``read`` and ``write`` methods onto registered classes.

        Will modify classes which have been registered to a reader or writer to have
        ``read`` and ``write`` methods which will contain documentation specifying
        useable formats for that class.

        The actual functionality will be a pass-through to :func:`skbio.io.read` and
        :func:`skbio.io.write` respectively.
        """
        reads = set()
        writes = set()
        for lookup in self._lookups:
            for format in lookup.values():
                reads |= format.monkey_patched_readers
                writes |= format.monkey_patched_writers

        for cls in reads:
            self._apply_read(cls)

        for cls in writes:
            self._apply_write(cls)

    def _apply_read(registry, cls):
        """Add read method if any formats have a reader for `cls`."""
        read_formats = registry.list_read_formats(cls)

        @classonlymethod
        def read(cls, file, format=None, **kwargs):
            return registry.read(file, into=cls, format=format, **kwargs)

        imports = registry._import_paths(read_formats)
        doc_list = registry._formats_for_docs(read_formats, imports)
        read.__func__.__doc__ = _read_docstring % {
            "name": cls.__name__,
            "list": doc_list,
            "see": "\n".join(imports),
        }
        cls.read = read

    def _apply_write(registry, cls):
        """Add write method if any formats have a writer for `cls`."""
        write_formats = registry.list_write_formats(cls)
        if not hasattr(cls, "default_write_format"):
            raise NotImplementedError(
                "Classes with registered writers must provide a "
                "`default_write_format`. Please add `default_write_format`"
                " to '%s'." % cls.__name__
            )

        def write(self, file, format=cls.default_write_format, **kwargs):
            return registry.write(self, into=file, format=format, **kwargs)

        imports = registry._import_paths(write_formats)
        doc_list = registry._formats_for_docs(write_formats, imports)
        write.__doc__ = _write_docstring % {
            "name": cls.__name__,
            "list": doc_list,
            "see": "\n".join(imports),
            "default": cls.default_write_format,
        }

        cls.write = write

    def _import_paths(self, formats):
        lines = []
        for fmt in formats:
            lines.append("skbio.io.format." + fmt)
        return lines

    def _formats_for_docs(self, formats, imports):
        lines = []
        for fmt, imp in zip(formats, imports):
            lines.append("- ``'%s'`` (:mod:`%s`)" % (fmt, imp))
        return "\n".join(lines)


_read_docstring = """Create a new ``%(name)s`` instance from a file.

This is a convenience method for :func:`skbio.io.registry.read`. For
more information about the I/O system in scikit-bio, please see
:mod:`skbio.io`.

Supported file formats include:

%(list)s

Parameters
----------
file : openable (filepath, URL, filehandle, etc.)
    The location to read the given `format`. Something that is
    understood by :func:`skbio.io.util.open`. Filehandles are not
    automatically closed, it is the responsibility of the caller.
format : str, optional
    The format must be a format name with a reader for ``%(name)s``.
    If a `format` is not provided or is None, it will attempt to
    guess the format.
kwargs : dict, optional
    Keyword arguments passed to :func:`skbio.io.registry.read` and
    the file format reader for ``%(name)s``.

Returns
-------
%(name)s
    A new instance.

See Also
--------
write
skbio.io.registry.read
skbio.io.util.open
%(see)s

"""

_write_docstring = """Write an instance of ``%(name)s`` to a file.

This is a convenience method for :func:`skbio.io.registry.write`.
For more information about the I/O system in scikit-bio, please
see :mod:`skbio.io`.

Supported file formats include:

%(list)s

Parameters
----------
file : openable (filepath, URL, filehandle, etc.)
    The location to write the given `format` into.  Something
    that is understood by :func:`skbio.io.util.open`. Filehandles
    are not automatically closed, it is the responsibility of the
    caller.
format : str
    The format must be a registered format name with a writer for
    ``%(name)s``.
    Default is `'%(default)s'`.
kwargs : dict, optional
    Keyword arguments passed to :func:`skbio.io.registry.write`
    and the file format writer.

See Also
--------
read
skbio.io.registry.write
skbio.io.util.open
%(see)s

"""


class Format:
    """Defines a format on which readers/writers/sniffer can be registered.

    Parameters
    ----------
    name : str
        The name of this format.
    encoding : str, optional
        What the default encoding of this format is. If set to 'binary' then
        all registered handlers will receive an :class:`io.BufferedReader` or
        :class:`io.BufferedWriter` instead of an :class:`io.TextIOBase`. The
        user will also be unable to override the encoding in that case.
    newline : str, optional
        What the default newline handling of this format is. Default is to use
        universal newline handling.

    """

    @property
    def name(self):
        """The name of this format."""
        return self._name

    @property
    def is_binary_format(self):
        """Return True if this is a binary format."""
        return self._encoding == "binary"

    @property
    def sniffer_function(self):
        """The sniffer function associated with this format."""
        return self._sniffer_function

    @property
    def readers(self):
        """Dictionary that maps classes to their readers for this format."""
        return self._readers

    @property
    def writers(self):
        """Dictionary that maps classes to their writers for this format."""
        return self._writers

    @property
    def monkey_patched_readers(self):
        """Set of classes bound to readers to monkey patch."""
        return self._monkey_patch["read"]

    @property
    def monkey_patched_writers(self):
        """Set of classes bound to writers to monkey patch."""
        return self._monkey_patch["write"]

    def __init__(self, name, encoding=None, newline=None):
        """Initialize format for registering sniffers, readers, and writers."""
        self._encoding = encoding
        self._newline = newline
        self._name = name

        self._sniffer_function = None
        self._readers = {}
        self._writers = {}
        self._monkey_patch = {"read": set(), "write": set()}

    def sniffer(self, override=False):
        r"""Decorate a function to act as the sniffer for this format.

        The function should take one argument which will be an implementation
        of either :class:`io.TextIOBase` or :class:`io.BufferedReader`
        depending on if the format is text or binary, respectively.

        The sniffer will always receive a filehandle which is pointing to the
        beginning of the file. It must return a tuple of bool and a dict of
        suggested keyword arguments (if any) to pass to the reader.

        .. note:: Keyword arguments are not permitted in `sniffers`.
           `Sniffers` may not raise exceptions; if an exception is thrown by a
           `sniffer`, the user will be asked to report it on our `issue tracker
           <https://github.com/scikit-bio/scikit-bio/issues/>`_.


        Parameters
        ----------
        override : bool, optional
            If True, the existing sniffer will be overriden.

        Raises
        ------
        DuplicateRegistrationError
            When `override` is False and a sniffer is already registered for
            this format.

        Examples
        --------
        >>> from skbio.io.registry import Format
        >>> # If developing a new format for skbio, use the create_format()
        >>> # factory instead of this constructor.
        >>> myformat = Format('myformat')
        >>> @myformat.sniffer()
        ... def myformat_sniffer(fh):
        ...     check = fh.read(8) == "myformat"
        ...     if check:
        ...         version = int(fh.read(1))
        ...         return True, {'version': version}
        ...     return False, {}
        ...
        >>> myformat_sniffer(["myformat2\n", "some content\n"])
        (True, {'version': 2})
        >>> myformat_sniffer(["something else\n"])
        (False, {})

        """
        if isinstance(override, bool) is not True:
            raise InvalidRegistrationError(
                "`override` must be a bool not %r" % override
            )

        if not override and self._sniffer_function is not None:
            raise DuplicateRegistrationError(
                "A sniffer is already registered" " to format: %s" % self._name
            )

        def decorator(sniffer):
            @wraps(sniffer)
            def wrapped_sniffer(
                file,
                encoding=self._encoding,
                errors="ignore",
                newline=self._newline,
                **kwargs,
            ):
                self._validate_encoding(encoding)
                if encoding == "binary":
                    # Errors is irrelevant so set to default to prevent raising
                    # a usage exception in open.
                    errors = _open_kwargs["errors"]
                with open_file(
                    file,
                    mode="r",
                    encoding=encoding,
                    newline=newline,
                    errors=errors,
                    **kwargs,
                ) as fh:
                    try:
                        # Some formats may have headers which indicate their
                        # format sniffers should be able to rely on the
                        # filehandle to point at the beginning of the file.
                        fh.seek(0)
                        return sniffer(fh)
                    except UnicodeDecodeError:
                        pass
                    except Exception:
                        warn(
                            f"'{sniffer.__name__}' has encountered a "
                            "problem.\nPlease send the following to our "
                            "issue tracker at\nhttps://github.com/scikit-"
                            "bio/scikit-bio/issues\n\n"
                            f"{traceback.format_exc()}",
                            FormatIdentificationWarning,
                        )

                    return False, {}

            self._sniffer_function = wrapped_sniffer
            return wrapped_sniffer

        return decorator

    def reader(self, cls, monkey_patch=True, override=False):
        r"""Decorate a function to act as the reader for a class in this format.

        The function should take an argument which will be an implementation
        of either :class:`io.TextIOBase` or :class:`io.BufferedReader`
        depending on if the format is text or binary, respectively. Any kwargs
        given by the user which are not handled by :func:`skbio.io.util.open`
        will be passed into the function. Any kwarg with a default of
        `FileSentinel` will transform user input for that parameter into a
        filehandle or `None` if not provided.

        Parameters
        ----------
        cls : type or None
            The class which the function will be registered to handle. If
            None, it is assumed that the function will produce a generator.
        monkey_patch : bool, optional
            Whether to allow an IORegistry to attach a `read` method to `cls`
            with this format listed as an option.
        override : bool, optional
            If True, any existing readers for `cls` in this format will be
            overriden.

        Raises
        ------
        DuplicateRegistrationError
            When `override` is False and a reader is already registered to
            `cls` for this format.

        Examples
        --------
        >>> from skbio.io.registry import Format, IORegistry
        >>> registry = IORegistry()
        >>> myformat = Format('myformat')
        >>> registry.add_format(myformat)
        >>> # If developing a new format for skbio, use the create_format()
        >>> # factory instead of the above.
        >>> class MyObject:
        ...     def __init__(self, content):
        ...         self.content = content
        ...
        >>> @myformat.reader(MyObject)
        ... def myformat_reader(fh):
        ...     return MyObject(fh.readlines()[1:])
        ...
        >>> registry.monkey_patch() # If developing skbio, this isn't needed
        >>> MyObject.read(["myformat2\n", "some content here!\n"],
        ...               format='myformat').content
        ['some content here!\n']

        """
        self._check_registration(cls)

        def decorator(reader_function):
            file_params = find_sentinels(reader_function, FileSentinel)
            # This split has to occur for the same reason as in IORegistry.read
            if cls is not None:

                @wraps(reader_function)
                def wrapped_reader(
                    file, encoding=self._encoding, newline=self._newline, **kwargs
                ):
                    file_keys, files, io_kwargs = self._setup_locals(
                        file_params, file, encoding, newline, kwargs
                    )
                    with open_files(files, mode="r", **io_kwargs) as fhs:
                        # The primary file is at the end of fh because append
                        # is cheaper than insert
                        kwargs.update(zip(file_keys, fhs[:-1]))
                        return reader_function(fhs[-1], **kwargs)
            else:

                @wraps(reader_function)
                def wrapped_reader(
                    file, encoding=self._encoding, newline=self._newline, **kwargs
                ):
                    file_keys, files, io_kwargs = self._setup_locals(
                        file_params, file, encoding, newline, kwargs
                    )
                    with open_files(files, mode="r", **io_kwargs) as fhs:
                        kwargs.update(zip(file_keys, fhs[:-1]))
                        yield from reader_function(fhs[-1], **kwargs)

            self._add_reader(cls, wrapped_reader, monkey_patch, override)
            return wrapped_reader

        return decorator

    def writer(self, cls, monkey_patch=True, override=False):
        r"""Decorate a function to act as the writer for a class in this format.

        The function should take an instance of `cls` as its first argument
        and the second argument is a filehandle which will be an implementation
        of either :class:`io.TextIOBase` or :class:`io.BufferedWriter`
        depending on if the format is text or binary, respectively. Any kwargs
        given by the user which are not handled by :func:`skbio.io.util.open`
        will be passed into the function. Any kwarg with a default of
        `FileSentinel` will transform user input for that parameter into a
        filehandle or `None` if not provided.

        Parameters
        ----------
        cls : type or None
            The class which the function will be registered to handle. If
            None, it is assumed that the function will consume a generator.
        monkey_patch : bool, optional
            Whether to allow an IORegistry to attach a `write` method to `cls`
            with this format listed as an option.
        override : bool, optional
            If True, any existing writers for `cls` in this format will be
            overriden.

        Raises
        ------
        DuplicateRegistrationError
            When `override` is False and a writer is already registered to
            `cls` for this format.

        Examples
        --------
        >>> from skbio.io.registry import Format, IORegistry
        >>> registry = IORegistry()
        >>> myformat = Format('myformat')
        >>> registry.add_format(myformat)
        >>> # If developing a new format for skbio, use the create_format()
        >>> # factory instead of the above.
        >>> class MyObject:
        ...     default_write_format = 'myformat'
        ...     def __init__(self, content):
        ...         self.content = content
        ...
        >>> @myformat.writer(MyObject)
        ... def myformat_reader(obj, fh):
        ...     fh.write("myformat2\n")
        ...     for c in obj.content:
        ...         fh.write(c)
        ...
        >>> registry.monkey_patch() # If developing skbio, this isn't needed
        >>> obj = MyObject(["some content here!\n"])
        >>> obj.write([], format='myformat')
        ['myformat2\n', 'some content here!\n']

        """
        self._check_registration(cls)

        def decorator(writer_function):
            file_params = find_sentinels(writer_function, FileSentinel)

            @wraps(writer_function)
            def wrapped_writer(
                obj, file, encoding=self._encoding, newline=self._newline, **kwargs
            ):
                file_keys, files, io_kwargs = self._setup_locals(
                    file_params, file, encoding, newline, kwargs
                )
                with open_files(files, mode="w", **io_kwargs) as fhs:
                    kwargs.update(zip(file_keys, fhs[:-1]))
                    writer_function(obj, fhs[-1], **kwargs)

            self._add_writer(cls, wrapped_writer, monkey_patch, override)
            return wrapped_writer

        return decorator

    def _check_registration(self, cls):
        if cls is not None and not inspect.isclass(cls):
            raise InvalidRegistrationError(
                "`cls` must be a class or None, not" " %r" % cls
            )

    def _setup_locals(self, file_params, file, encoding, newline, kwargs):
        self._validate_encoding(encoding)
        io_kwargs = self._pop_io_kwargs(kwargs, encoding, newline)
        file_keys, files = self._setup_file_args(kwargs, file_params)
        files.append(file)

        return file_keys, files, io_kwargs

    def _validate_encoding(self, encoding):
        if encoding != self._encoding:
            if self._encoding == "binary":
                raise ValueError("Encoding must be 'binary' for %r" % self.name)
            if encoding == "binary":
                raise ValueError("Encoding must not be 'binary' for %r" % self.name)

    def _pop_io_kwargs(self, kwargs, encoding, newline):
        io_kwargs = dict(encoding=encoding, newline=newline)
        for key in _open_kwargs:
            if key in kwargs:
                io_kwargs[key] = kwargs.pop(key)
        return io_kwargs

    def _setup_file_args(self, kwargs, file_params):
        file_keys = []
        files = []
        for param in file_params:
            arg = kwargs.get(param, None)
            if arg is not None:
                file_keys.append(param)
                files.append(arg)
            else:
                # set to None to mask FileSentinel when user neglected argument
                kwargs[param] = None

        return file_keys, files

    def _add_writer(self, cls, writer, monkey_patch, override):
        if cls in self._writers and not override:
            raise DuplicateRegistrationError(
                "There is already a writer"
                " registered to %s in format: %s" % (cls, self._name)
            )
        self._writers[cls] = writer
        if monkey_patch and cls is not None:
            self._monkey_patch["write"].add(cls)

    def _add_reader(self, cls, reader, monkey_patch, override):
        if cls in self._readers and not override:
            raise DuplicateRegistrationError(
                "There is already a reader"
                " registered to %s in format: %s" % (cls, self._name)
            )
        self._readers[cls] = reader
        if monkey_patch and cls is not None:
            self._monkey_patch["read"].add(cls)


io_registry = IORegistry()


@wraps(IORegistry.sniff)
def sniff(file, **kwargs):
    """Detect the format of a given `file` and suggest kwargs."""
    return io_registry.sniff(file, **kwargs)


@wraps(IORegistry.read)
def read(file, format=None, into=None, verify=True, **kwargs):
    """Read data from file using specified format."""
    return io_registry.read(file, format=format, into=into, verify=verify, **kwargs)


@wraps(IORegistry.write)
def write(obj, format, into, **kwargs):
    """Write data to file using specified format."""
    return io_registry.write(obj, format, into, **kwargs)


@wraps(IORegistry.create_format)
def create_format(*args, **kwargs):
    """Make a new format."""
    return io_registry.create_format(*args, **kwargs)
