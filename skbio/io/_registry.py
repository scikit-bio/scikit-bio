# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from warnings import warn
import types
import traceback
import inspect

from future.builtins import zip

from . import (UnrecognizedFormatError, InvalidRegistrationError,
               DuplicateRegistrationError, ArgumentOverrideWarning,
               FormatIdentificationWarning)
from .util import open_files, reopen_file

_formats = {}
_sniffers = {}
_aliases = {}
_empty_file_format = '<emptyfile>'

# We create a class and instantiate it dynamically so that exceptions are more
# obvious and so that only one object exists without copying this line.
FileSentinel = type('FileSentinel', (object, ), {})()


def _override_kwargs(kw, fmt_kw, warn_user):
    for key in kw:
        if key in fmt_kw and fmt_kw[key] != kw[key] and warn_user:
            warn('Best guess was: %s=%s, continuing with user supplied: %s' % (
                key, str(fmt_kw[key]), str(kw[key])
            ), ArgumentOverrideWarning)
        fmt_kw[key] = kw[key]
    return fmt_kw


def register_sniffer(format, binary=False):
    """Return a decorator for a sniffer function.

    A decorator factory for sniffer functions. Sniffers may only be registered
    to simple formats. Sniffers for compound formats are automatically
    generated from their component simple formats.

    A sniffer function should have at least the following signature:
    ``<format_name>_sniffer(fh)``. `fh` is **always** an open filehandle.
    This decorator provides the ability to use filepaths in the same argument
    position as `fh`. They will automatically be opened and closed.

    **The sniffer must not close the filehandle**, cleanup will be
    handled external to the sniffer and is not its concern.

    `**kwargs` are not passed to a sniffer, and a sniffer must not use them.

    The job of a sniffer is to determine if a file appears to be in the given
    format and to 'sniff' out any kwargs that would be of use to a reader
    function.

    The sniffer **must** return a tuple of (True, <kwargs dict>) if it believes
    `fh` is a given `format`. Otherwise it should return (False, {}).

    .. note:: Failure to adhere to the above interface specified for a sniffer
       will result in unintended side-effects.

    The sniffer may determine membership of a file in as many or as few
    lines of the file as it deems necessary.

    Parameters
    ----------
    format : str
        A format name which a decorated sniffer will be bound to.

    Returns
    -------
    function
        A decorator to be used on a sniffer. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists a
        *sniffer* bound to the `format`.

    See Also
    --------
    skbio.io.sniff

    """
    def decorator(sniffer):
        if format in _sniffers:
            raise DuplicateRegistrationError(
                msg="'%s' already has a sniffer." % format)

        def wrapped_sniffer(fp, mode='r', binary=binary, **kwargs):
            with reopen_file(fp, binary=binary) as fh:
                try:
                    return sniffer(fh, **kwargs)
                except Exception:
                    warn("'%s' has encountered a problem.\n"
                         "Please send the following to our issue tracker at\n"
                         "https://github.com/biocore/scikit-bio/issues\n\n"
                         "%s" % (sniffer.__name__, traceback.format_exc()),
                         FormatIdentificationWarning)
                    return False, {}

        wrapped_sniffer.__doc__ = sniffer.__doc__
        wrapped_sniffer.__name__ = sniffer.__name__

        _sniffers[format] = wrapped_sniffer
        return wrapped_sniffer
    return decorator


def register_reader(format, cls=None, binary=False):
    """Return a decorator for a reader function.

    A decorator factory for reader functions.

    A reader function should have at least the following signature:
    ``<format_name>_to_<class_name_or_generator>(fh)``. `fh` is **always** an
    open filehandle. This decorator provides the ability to use filepaths in
    the same argument position as `fh`. They will automatically be opened and
    closed.

    **The reader must not close the filehandle**, cleanup will be
    handled external to the reader and is not its concern. This is true even
    in the case of generators.

    Any additional `**kwargs` will be passed to the reader and may
    be used if necessary.

    The reader **must** return an instance of `cls` if `cls` is not None.
    Otherwise the reader must return a generator. The generator need not deal
    with closing the `fh`. That is already handled by this decorator.

    .. note:: Failure to adhere to the above interface specified for a reader
       will result in unintended side-effects.

    Parameters
    ----------
    format : str
        A format name which a decorated reader will be bound to.
    cls : type, optional
        The class which a decorated reader will be bound to. When `cls` is None
        the reader will be bound as returning a generator.
        Default is None.

    Returns
    -------
    function
        A decorator to be used on a reader. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists a
        *reader* bound to the same permutation of `fmt` and `cls`.

    See Also
    --------
    skbio.io.read

    """
    def decorator(reader):
        format_class = _formats.setdefault(format, {}).setdefault(cls, {})

        if 'reader' in format_class:
            raise DuplicateRegistrationError('reader', format, cls)

        file_args = []
        reader_spec = inspect.getargspec(reader)
        if reader_spec.defaults is not None:
            # Concept from http://stackoverflow.com/a/12627202/579416
            for key, default in zip(
                    reader_spec.args[-len(reader_spec.defaults):],
                    reader_spec.defaults):
                if default is FileSentinel:
                    file_args.append(key)

        # We wrap the reader so that basic file handling can be managed
        # externally from the business logic.
        if cls is None:
            def wrapped_reader(fp, mode='r', binary=binary,
                               mutate_fh=False, **kwargs):
                file_keys = []
                files = [fp]
                for file_arg in file_args:
                    if file_arg in kwargs:
                        if kwargs[file_arg] is not None:
                            file_keys.append(file_arg)
                            files.append(kwargs[file_arg])
                    else:
                        kwargs[file_arg] = None

                # Do we want binary argument to come from the user
                # (wrapped_reader) or from the reader (register_reader)?
                with open_files(files, mode=mode, binary=binary) as fhs:
                    try:
                        for key, fh in zip(file_keys, fhs[1:]):
                            kwargs[key] = fh

                        generator = reader(fhs[0], **kwargs)
                        if not isinstance(generator, types.GeneratorType):
                            # Raise an exception to be handled next line,
                            # because although reader executed without error,
                            # it is not a generator.
                            raise Exception()
                    # If an exception is thrown at this point, it cannot
                    # be a generator. If there was a `yield` statment, then
                    # Python would have returned a generator regardless of the
                    # content. This does not preclude the generator from
                    # throwing exceptions.
                    except Exception:
                            raise InvalidRegistrationError("'%s' is not a "
                                                           "generator." %
                                                           reader.__name__)

                    while True:
                        yield next(generator)

        else:
            # When an object is instantiated we don't need to worry about the
            # original position at every step, only at the end.
            def wrapped_reader(fp, mode='r', binary=binary,
                               mutate_fh=False, **kwargs):
                file_keys = []
                files = [fp]
                for file_arg in file_args:
                    if file_arg in kwargs:
                        if kwargs[file_arg] is not None:
                            file_keys.append(file_arg)
                            files.append(kwargs[file_arg])
                    else:
                        kwargs[file_arg] = None

                with open_files(files, mode=mode, binary=binary) as fhs:
                    for key, fh in zip(file_keys, fhs[1:]):
                        kwargs[key] = fh
                    return reader(fhs[0], **kwargs)

        wrapped_reader.__doc__ = reader.__doc__
        wrapped_reader.__name__ = reader.__name__

        format_class['reader'] = wrapped_reader
        return wrapped_reader
    return decorator


def register_writer(format, cls=None, binary=False):
    """Return a decorator for a writer function.

    A decorator factory for writer functions.

    A writer function should have at least the following signature:
    ``<class_name_or_generator>_to_<format_name>(obj, fh)``. `fh` is **always**
    an open filehandle. This decorator provides the ability to use filepaths in
    the same argument position as `fh`. They will automatically be opened and
    closed.

    **The writer must not close the filehandle**, cleanup will be
    handled external to the reader and is not its concern.

    Any additional `**kwargs` will be passed to the writer and may be used if
    necessary.

    The writer must not return a value. Instead it should only mutate the `fh`
    in a way consistent with it's purpose.

    If the writer accepts a generator, it should exhaust the generator to
    ensure that the potentially open filehandle backing said generator is
    closed.

    .. note:: Failure to adhere to the above interface specified for a writer
       will result in unintended side-effects.

    Parameters
    ----------
    format : str
        A format name which a decorated writer will be bound to.
    cls : type, optional
        The class which a decorated writer will be bound to. If `cls` is None
        the writer will be bound as expecting a generator.
        Default is None.

    Returns
    -------
    function
        A decorator to be used on a writer. The decorator will raise a
        ``skbio.io.DuplicateRegistrationError`` if there already exists a
        *writer* bound to the same permutation of `fmt` and `cls`.

    See Also
    --------
    skbio.io.write
    skbio.io.get_writer

    """
    def decorator(writer):
        format_class = _formats.setdefault(format, {}).setdefault(cls, {})

        if 'writer' in format_class:
            raise DuplicateRegistrationError('writer', format, cls)

        file_args = []
        writer_spec = inspect.getargspec(writer)
        if writer_spec.defaults is not None:
            # Concept from http://stackoverflow.com/a/12627202/579416
            for key, default in zip(
                    writer_spec.args[-len(writer_spec.defaults):],
                    writer_spec.defaults):
                if default is FileSentinel:
                    file_args.append(key)

        # We wrap the writer so that basic file handling can be managed
        # externally from the business logic.
        def wrapped_writer(obj, fp, mode='w', binary=binary,
                           gzip=False, compresslevel=9, **kwargs):
            file_keys = []
            files = [fp]
            for file_arg in file_args:
                if file_arg in kwargs:
                    if kwargs[file_arg] is not None:
                        file_keys.append(file_arg)
                        files.append(kwargs[file_arg])
                else:
                    kwargs[file_arg] = None

            with open_files(files, mode=mode, binary=binary,
                            gzip=gzip, compresslevel=compresslevel) as fhs:
                for key, fh in zip(file_keys, fhs[1:]):
                    kwargs[key] = fh
                writer(obj, fhs[0], **kwargs)

        wrapped_writer.__doc__ = writer.__doc__
        wrapped_writer.__name__ = writer.__name__

        format_class['writer'] = wrapped_writer
        return wrapped_writer
    return decorator


def list_read_formats(cls):
    """Return a list of available read formats for a given `cls` type.

    Parameters
    ----------
    cls : type
        The class which will be used to determine what read formats exist for
        an instance of `cls`.

    Returns
    -------
    list
        A list of available read formats for an instance of `cls`. List may be
        empty.

    See Also
    --------
    skbio.io.register_reader

    """
    return _rw_list_formats('reader', cls)


def list_write_formats(cls):
    """Return a list of available write formats for a given `cls` instance.

    Parameters
    ----------
    cls : type
        The class which will be used to determine what write formats exist for
        an instance of `cls`.

    Returns
    -------
    list
        A list of available write formats for an instance of `cls`. List may be
        empty.

    See Also
    --------
    skbio.io.register_writer

    """
    return _rw_list_formats('writer', cls)


def _rw_list_formats(name, cls):
    formats = []
    for fmt in _formats:
        if cls in _formats[fmt] and name in _formats[fmt][cls]:
            formats.append(fmt)
    return formats


def get_sniffer(format):
    """Return a sniffer for a format.

    Parameters
    ----------
    format : str
        A format string which has a registered sniffer.

    Returns
    -------
    function or None
        Returns a sniffer function if one exists for the given `fmt`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_sniffer

    """
    return _sniffers.get(format, None)


def get_reader(format, cls=None):
    """Return a reader for a format.

    Parameters
    ----------
    format : str
        A registered format string.
    cls : type, optional
        The class which the reader will return an instance of. If `cls` is
        None, the reader will return a generator.
        Default is None.

    Returns
    -------
    function or None
        Returns a reader function if one exists for a given `fmt` and `cls`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_reader

    """
    return _rw_getter('reader', format, cls)


def get_writer(format, cls=None):
    """Return a writer for a format.

    Parameters
    ----------
    format : str
        A registered format string.
    cls : type, optional
        The class which the writer will expect an instance of. If `cls` is
        None, the writer will expect a generator that is identical to what
        is returned by ``get_reader(<some_format>, None)``.
        Default is None.

    Returns
    -------
    function or None
        Returns a writer function if one exists for a given `fmt` and `cls`.
        Otherwise it will return None.

    See Also
    --------
    skbio.io.register_writer
    skbio.io.get_reader

    """
    return _rw_getter('writer', format, cls)


def _rw_getter(name, fmt, cls):
    if fmt in _formats:
        if cls in _formats[fmt] and name in _formats[fmt][cls]:
                return _formats[fmt][cls][name]
    return None


def sniff(fp, cls=None, mode='r'):
    """Attempt to guess the format of a file and return format str and kwargs.

    Parameters
    ----------
    fp : filepath or filehandle
        The provided file to guess the format of. Filepaths and HTTP/HTTPS URLs
        are automatically closed; filehandles are the responsibility of the
        caller.
    cls : type, optional
        A provided class that restricts the search for the format. Only formats
        which have a registered reader or writer for the given `cls` will be
        tested.
        Default is None.

    Returns
    -------
    (str, kwargs)
        A format name and kwargs for the corresponding reader.

    Raises
    ------
    UnrecognizedFormatError
        This occurs when the format is not 'claimed' by any registered sniffer
        or when the format is ambiguous and has been 'claimed' by more than one
        sniffer.

    See Also
    --------
    skbio.io.register_sniffer
    skbio.io.util.open_file

    """
    possibles = []
    for fmt in _sniffers:
        if cls is not None and fmt != _empty_file_format and (
                fmt not in _formats or cls not in _formats[fmt]):
            continue
        format_sniffer = _sniffers[fmt]
        is_format, fmt_kwargs = format_sniffer(fp, mode=mode)
        if is_format:
            possibles.append(fmt)
            kwargs = fmt_kwargs

    if not possibles:
        raise UnrecognizedFormatError("Cannot guess the format for %s."
                                      % str(fp))
    if len(possibles) > 1:
        raise UnrecognizedFormatError("File format is ambiguous, may be"
                                      " one of %s." % str(possibles))
    return possibles[0], kwargs


def read(fp, format=None, into=None, verify=True, mode='r', **kwargs):
    """Read a supported skbio file format into an instance or a generator.

    This function is able to reference and execute all *registered* read
    operations in skbio.

    Parameters
    ----------
    fp : filepath or filehandle
        The location to read the given `format` `into`. Filepaths or HTTP/HTTPS
        URLs are automatically closed when read; filehandles are the
        responsibility of the caller. In the case of a generator, a filepath
        will be closed when ``StopIteration`` is raised; filehandles are still
        the responsibility of the caller.
    format : str, optional
        The format must be a format name with a reader for the given
        `into` class. If a `format` is not provided or is None, all
        registered sniffers for the provied `into` class will be evaluated to
        attempt to guess the format.
        Default is None.
    into : type, optional
        A class which has a registered reader for a given `format`. If `into`
        is not provided or is None, read will return a generator.
        Default is None.
    verify : bool, optional
        Whether or not to confirm the format of a file if `format` is provided.
        Will raise a ``skbio.io.FormatIdentificationWarning`` if the sniffer of
        `format` returns False.
        Default is True.
    mode : str, optional
        The read mode. This is passed to `open(fp, mode)` internally.
        Default is 'U'
    kwargs : dict, optional
        Will be passed directly to the appropriate reader.

    Returns
    -------
    object or generator
        If `into` is not None, an instance of the `into` class will be
        provided with internal state consistent with the provided file.
        If `into` is None, a generator will be returned.

    Raises
    ------
    ValueError
        Raised when `format` and `into` are both None.
    skbio.io.UnrecognizedFormatError
        Raised when a reader could not be found for a given `format` or the
        format could not be guessed.
    skbio.io.FormatIdentificationWarning
        Raised when `verify` is True and the sniffer of a `format` provided a
        kwarg value that did not match the user's kwarg value.

    See Also
    --------
    skbio.io.register_reader
    skbio.io.register_sniffer
    skbio.io.util.open_files

    """
    if format is None and into is None:
        raise ValueError("`format` and `into` cannot both be None.")

    if format is None:
        format, fmt_kwargs = sniff(fp, cls=into, mode=mode)
        kwargs = _override_kwargs(kwargs, fmt_kwargs, verify)
    elif verify:
        sniffer = get_sniffer(format)
        if sniffer is not None:
            is_format, fmt_kwargs = sniffer(fp)
            if not is_format:
                warn("%s could not be positively identified as %s file." %
                     (str(fp), format),
                     FormatIdentificationWarning)
            else:
                kwargs = _override_kwargs(kwargs, fmt_kwargs, True)

    reader = get_reader(format, into)
    if reader is None:
        raise UnrecognizedFormatError("Cannot read %s into %s, no reader "
                                      "found." % (format, into.__name__
                                                  if into is not None
                                                  else 'generator'))
    return reader(fp, mode=mode, **kwargs)


def write(obj, format, into, mode='w', gzip=False, compresslevel=9, **kwargs):
    """Write a supported skbio file format from an instance or a generator.

    This function is able to reference and execute all *registered* write
    operations in skbio.

    Parameters
    ----------
    obj : object
        The object must have a registered writer for a provided `format`.
    format : str
        The format must be a registered format name with a writer for the given
        `obj`.
    into : filepath or filehandle
        The location to write the given `format` from `obj` into. Filepaths are
        automatically closed when written; filehandles are the responsibility
        of the caller.
    mode : str, optional
        The write mode. This is passed to `open(fp, mode)` internally.
        Default is 'w'.
    kwargs : dict, optional
        Will be passed directly to the appropriate writer.

    Raises
    ------
    skbio.io.UnrecognizedFormatError
        Raised when a writer could not be found for the given `format` and
        `obj`.

    See Also
    --------
    skbio.io.register_writer

    """
    cls = None
    if not isinstance(obj, types.GeneratorType):
        cls = obj.__class__
    writer = get_writer(format, cls)
    if writer is None:
        raise UnrecognizedFormatError("Cannot write %s into %s, no %s writer "
                                      "found." % (format, str(into),
                                                  'generator' if cls is None
                                                  else str(cls)))

    writer(obj, into, mode=mode, gzip=gzip,
           compresslevel=compresslevel, **kwargs)


# This is meant to be a handy indicator to the user that they have done
# something wrong.
@register_sniffer(_empty_file_format)
def empty_file_sniffer(fh):
    for line in fh:
        if line.strip():
            return False, {}
    return True, {}


def initialize_oop_interface():
    classes = set()
    # Find each potential class
    for fmt in _formats:
        for cls in _formats[fmt]:
            classes.add(cls)
    # Add readers and writers for each class
    for cls in classes:
        if cls is not None:
            _apply_read(cls)
            _apply_write(cls)


def _apply_read(cls):
    """Add read method if any formats have a registered reader for `cls`."""
    skbio_io_read = globals()['read']
    read_formats = list_read_formats(cls)
    if read_formats:
        @classmethod
        def read(cls, fp, format=None, **kwargs):
            return skbio_io_read(fp, into=cls, format=format, **kwargs)

        read.__func__.__doc__ = _read_docstring % (
            cls.__name__,
            _formats_for_docs(read_formats),
            cls.__name__,
            cls.__name__,
            cls.__name__,
            _import_paths(read_formats)
        )
        cls.read = read


def _apply_write(cls):
    """Add write method if any formats have a registered writer for `cls`."""
    skbio_io_write = globals()['write']
    write_formats = list_write_formats(cls)
    if write_formats:
        if not hasattr(cls, 'default_write_format'):
            raise NotImplementedError(
                "Classes with registered writers must provide a "
                "`default_write_format`. Please add `default_write_format` to"
                " '%s'." % cls.__name__)

        def write(self, fp, format=cls.default_write_format, **kwargs):
            skbio_io_write(self, into=fp, format=format, **kwargs)

        write.__doc__ = _write_docstring % (
            cls.__name__,
            _formats_for_docs(write_formats),
            cls.__name__,
            cls.default_write_format,
            _import_paths(write_formats)
        )
        cls.write = write


def _import_paths(formats):
    lines = []
    for fmt in formats:
        lines.append("skbio.io." + fmt)
    return '\n'.join(lines)


def _formats_for_docs(formats):
    lines = []
    for fmt in formats:
        lines.append("- ``'%s'`` (:mod:`skbio.io.%s`)" % (fmt, fmt))
    return '\n'.join(lines)


_read_docstring = """Create a new ``%s`` instance from a file.

This is a convenience method for :mod:`skbio.io.read`. For more
information about the I/O system in scikit-bio, please see
:mod:`skbio.io`.

Supported file formats include:

%s

Parameters
----------
fp : filepath or filehandle
    The location to read the given `format`. Filepaths are
    automatically closed when read; filehandles are the
    responsibility of the caller.
format : str, optional
    The format must be a format name with a reader for ``%s``.
    If a `format` is not provided or is None, it will attempt to
    guess the format.
kwargs : dict, optional
    Keyword arguments passed to :mod:`skbio.io.read` and the file
    format reader for ``%s``.

Returns
-------
%s
    A new instance.

See Also
--------
write
skbio.io.read
%s

"""

_write_docstring = """Write an instance of ``%s`` to a file.

This is a convenience method for :mod:`skbio.io.write`. For more
information about the I/O system in scikit-bio, please see
:mod:`skbio.io`.

Supported file formats include:

%s

Parameters
----------
fp : filepath or filehandle
    The location to write the given `format` into. Filepaths are
    automatically closed when written; filehandles are the
    responsibility of the caller.
format : str
    The format must be a registered format name with a writer for
    ``%s``.
    Default is `'%s'`.
kwargs : dict, optional
    Keyword arguments passed to :mod:`skbio.io.write` and the
    file format writer.

See Also
--------
read
skbio.io.write
%s

"""
