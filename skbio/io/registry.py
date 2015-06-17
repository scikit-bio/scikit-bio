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
import copy
import traceback
from functools import wraps

from future.builtins import zip

from . import (UnrecognizedFormatError, InvalidRegistrationError,
               DuplicateRegistrationError, ArgumentOverrideWarning,
               FormatIdentificationWarning)
from .util import resolve_file, open_file, open_files, _d as _open_kwargs
from skbio.util._misc import Inherit, make_sentinel, find_sentinels
import skbio.io

FileSentinel = make_sentinel("FileSentinel")


class IORegistry(object):
    def __init__(self):
        self.binary_formats = {}
        self.text_formats = {}
        self._lookups = (self.binary_formats, self.text_formats)

    def add_format(self, format_object):
        if format_object.is_binary_format:
            self.binary_formats[format_object.name] = format_object
        else:
            self.text_formats[format_object.name] = format_object

    def get_sniffer(self, format_name):
        for lookup in self._lookups:
            if format_name in lookup:
                sniffer = lookup[format_name].sniffer_function
                if sniffer is not None:
                    return sniffer
        return None

    def get_reader(self, format_name, cls=None):
        return self._get_rw(format_name, cls, 'readers')

    def get_writer(self, format_name, cls=None):
        return self._get_rw(format_name, cls, 'writers')

    def _get_rw(self, format_name, cls, lookup_name):
        for lookup in self._lookups:
            if format_name in lookup:
                format_lookup = getattr(lookup[format_name], lookup_name)
                if cls in format_lookup:
                    return format_lookup[cls]
        return None

    def list_read_formats(self, cls):
        return list(self._iter_rw_formats(cls, 'readers'))

    def list_write_formats(self, cls):
        return list(self._iter_rw_formats(cls, 'writers'))

    def _iter_rw_formats(self, cls, lookup_name):
        for lookup in self._lookups:
            for format in lookup.values():
                if cls in getattr(format, lookup_name):
                    yield format.name

    def sniff(self, file, **kwargs):
        with resolve_file(file, mode='r', **kwargs) as (fh, is_binary_file):
            backup = fh.tell()
            if is_binary_file:
                matches = self._find_matches(fh, self.binary_formats, **kwargs)
            else:
                matches = []

            matches += self._find_matches(fh, self.text_formats, **kwargs)
            fh.seek(backup)

        if len(matches) > 1:
            raise UnrecognizedFormatError("File format for %r is ambiguous,"
                                          " may be one of: %r"
                                          % (file, [m for m, _ in matches]))
        elif len(matches) == 0:
            raise UnrecognizedFormatError("Could not detect the format of %r"
                                          % file)

        return matches[0]

    def _find_matches(self, file, lookup, **kwargs):
        matches = []
        for format in lookup.values():
            if format.sniffer_function is not None:
                file.seek(0)
                result = format.sniffer_function(file, **kwargs)

                if result[0]:
                    matches.append((format.name, result[1]))
        return matches

    def read(self, file, format=None, into=None, verify=True, **kwargs):
        io_kwargs = {k:kwargs[k] for k in _open_kwargs if k in kwargs}
        with resolve_file(file, **io_kwargs) as (file, _):
            if format is None:
                format, skwargs = self.sniff(file, **io_kwargs)
            else:
                if verify:
                    sniffer = self.get_sniffer(format)
                    if sniffer is not None:
                        backup = file.tell()
                        is_format, skwargs = sniffer(file, **io_kwargs)
                        file.seek(backup)
                        if not is_format:
                            warn("%r does not look like a %s file"
                                 % (file, format), FormatIdentificationWarning)

                    for key in skwargs:
                        if key not in kwargs:
                            kwargs[key] = skwargs[key]
                        elif kwargs[key] != skwargs[key]:
                            warn('Best guess was: %s=%r, continuing with user'
                                 ' supplied: %r' % (key, skwargs[key],
                                                    kwargs[key]),
                                 ArgumentOverrideWarning)

            reader = self.get_reader(format, into)
            if reader is None:
                raise UnrecognizedFormatError(
                    "Cannot read %s from %r, no %s reader found." %
                    (format, file, into.__name__ if into else 'generator'))

            return reader(file, **kwargs)

    def write(self, obj, format, into, **kwargs):
        cls = None
        if not isinstance(obj, types.GeneratorType):
            cls = obj.__class__
        writer = self.get_writer(format, cls)
        if writer is None:
            raise UnrecognizedFormatError(
                "Cannot write %r into %r, no %s writer found." %
                (format, into, obj.__class__.__name__))

        writer(obj, into, **kwargs)
        return into


    def monkey_patch(self, force=False):
        reads = set()
        writes = set()
        for lookup in self._lookups:
            for format in lookup.values():
                reads |= format.monkey_patch['read']
                writes |= format.monkey_patch['write']

        for cls in reads:
            self._apply_read(cls)

        for cls in writes:
            self._apply_write(cls)


    def _apply_read(registry, cls):
        """Add read method if any formats have a registered reader for `cls`."""
        read_formats = registry.list_read_formats(cls)
        if read_formats:
            @classmethod
            def read(cls, file, format=None, **kwargs):
                return registry.read(file, into=cls, format=format, **kwargs)

            imports = registry._import_paths(read_formats)
            doc_list = registry._formats_for_docs(read_formats, imports)
            read.__func__.__doc__ = _read_docstring % {
                'name': cls.__name__,
                'list': doc_list,
                'see': '\n'.join(imports)
            }
            cls.read = read


    def _apply_write(registry, cls):
        """Add write method if any formats have a registered writer for `cls`."""
        write_formats = registry.list_write_formats(cls)
        if write_formats:
            if not hasattr(cls, 'default_write_format'):
                raise NotImplementedError(
                    "Classes with registered writers must provide a "
                    "`default_write_format`. Please add `default_write_format` to"
                    " '%s'." % cls.__name__)

            def write(self, file, format=cls.default_write_format, **kwargs):
                return registry.write(self, into=file, format=format, **kwargs)


            imports = registry._import_paths(write_formats)
            doc_list = registry._formats_for_docs(write_formats, imports)
            write.__doc__ = _write_docstring % {
                'name': cls.__name__,
                'list': doc_list,
                'see': '\n'.join(imports),
                'default': cls.default_write_format
            }

            cls.write = write


    def _import_paths(self, formats):
        lines = []
        for fmt in formats:
            lines.append("skbio.io.formats." + fmt)
        return lines


    def _formats_for_docs(self, formats, imports):
        lines = []
        for fmt, imp in zip(formats, imports):
            lines.append("- ``'%s'`` (:mod:`%s`)" % (fmt, imp))
        return '\n'.join(lines)


_read_docstring = """Create a new ``%(name)s`` instance from a file.

This is a convenience method for :mod:`skbio.io.read`. For more
information about the I/O system in scikit-bio, please see
:mod:`skbio.io`.

Supported file formats include:

%(list)s

Parameters
----------
file : filepath or filehandle
    The location to read the given `format`. Filepaths are
    automatically closed when read; filehandles are the
    responsibility of the caller.
format : str, optional
    The format must be a format name with a reader for ``%(name)s``.
    If a `format` is not provided or is None, it will attempt to
    guess the format.
kwargs : dict, optional
    Keyword arguments passed to :mod:`skbio.io.read` and the file
    format reader for ``%(name)s``.

Returns
-------
%(name)s
    A new instance.

See Also
--------
write
skbio.io.read
%(see)s

"""

_write_docstring = """Write an instance of ``%(name)s`` to a file.

This is a convenience method for :mod:`skbio.io.write`. For more
information about the I/O system in scikit-bio, please see
:mod:`skbio.io`.

Supported file formats include:

%(list)s

Parameters
----------
file : filepath or filehandle
    The location to write the given `format` into. Filepaths are
    automatically closed when written; filehandles are the
    responsibility of the caller.
format : str
    The format must be a registered format name with a writer for
    ``%(name)s``.
    Default is `'%(default)s'`.
kwargs : dict, optional
    Keyword arguments passed to :mod:`skbio.io.write` and the
    file format writer.

See Also
--------
read
skbio.io.write
%(see)s

"""


class Format(object):
    def __init__(self, name, encoding=None, newline=None):
        self._encoding = encoding
        self._newline = newline
        self.name = name

        self.is_binary_format = encoding == 'binary'
        self.sniffer_function = None
        self.readers = {}
        self.writers = {}
        self.monkey_patch = {'read':set(), 'write': set()}

    def sniffer(self, sniffer):
        @wraps(sniffer)
        def wrapped_sniffer(file, encoding=self._encoding, errors='ignore',
                            newline=self._newline, **kwargs):
            self._validate_encoding(encoding)
            if encoding == 'binary':
                # Errors is irrelevant so set to default to prevent raising
                # a usage exception in open.
                errors = _open_kwargs['errors']
            with open_file(file, mode='r', encoding=encoding, newline=newline,
                           errors=errors, **kwargs) as fh:
                try:
                    return sniffer(fh)
                except UnicodeError:
                   pass # We don't need to warn in the event of a decoding
                        # error, but obviously we can't resolve it as True.
                except Exception:
                    warn("'%s' has encountered a problem.\n"
                         "Please send the following to our issue tracker at\n"
                         "https://github.com/biocore/scikit-bio/issues\n\n"
                         "%s" % (sniffer.__name__, traceback.format_exc()),
                         FormatIdentificationWarning)

                return False, {}

        self.sniffer_function = wrapped_sniffer
        return wrapped_sniffer

    def reader(self, cls=None, monkey_patch=True):
        def decorator(reader_function):
            file_params = find_sentinels(reader_function, FileSentinel)
            if cls is not None:
                @wraps(reader_function)
                def wrapped_reader(file, encoding=self._encoding,
                                   newline=self._newline, **kwargs):
                    file_keys, files, io_kwargs = self._setup_locals(
                        file_params, file, encoding, newline, kwargs)
                    with open_files(files, mode='r', **io_kwargs) as fhs:
                        kwargs.update(zip(file_keys, fhs[:-1]))
                        return reader_function(fhs[-1], **kwargs)
            else:
                @wraps(reader_function)
                def wrapped_reader(file, encoding=self._encoding,
                                   newline=self._newline, **kwargs):
                    file_keys, files, io_kwargs = self._setup_locals(
                        file_params, file, encoding, newline, kwargs)
                    with open_files(files, mode='r', **io_kwargs) as fhs:
                        kwargs.update(zip(file_keys, fhs[:-1]))
                        generator = reader_function(fhs[-1], **kwargs)
                        while True:
                            yield next(generator)

            self._add_reader(cls, wrapped_reader, monkey_patch)
            return wrapped_reader
        return decorator

    def writer(self, cls=None, monkey_patch=True):
        def decorator(writer_function):
            file_params = find_sentinels(writer_function, FileSentinel)

            @wraps(writer_function)
            def wrapped_writer(obj, file, encoding=self._encoding,
                               newline=self._newline, **kwargs):
                file_keys, files, io_kwargs = self._setup_locals(
                    file_params, file, encoding, newline, kwargs)
                with open_files(files, mode='w', **io_kwargs) as fhs:
                    kwargs.update(zip(file_keys, fhs[:-1]))
                    writer_function(obj, fhs[-1], **kwargs)

            self._add_writer(cls, wrapped_writer, monkey_patch)
            return wrapped_writer
        return decorator

    def _setup_locals(self, file_params, file, encoding, newline, kwargs):
        self._validate_encoding(encoding)
        io_kwargs = self._pop_io_kwargs(kwargs, encoding, newline)
        file_keys, files = self._setup_file_args(kwargs, file_params)
        files.append(file)

        return file_keys, files, io_kwargs

    def _validate_encoding(self, encoding):
        if encoding != self._encoding:
            if self._encoding == 'binary':
                raise ValueError("Encoding must be 'binary' for %r"
                                 % self.name)
            if encoding == 'binary':
                raise ValueError("Encoding must not be 'binary' for %r"
                                 % self.name)

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

    def _add_writer(self, cls, writer, monkey_patch):
        self.writers[cls] = writer
        if monkey_patch and cls is not None:
            self.monkey_patch['write'].add(cls)

    def _add_reader(self, cls, reader, monkey_patch):
        self.readers[cls] = reader
        if monkey_patch and cls is not None:
            self.monkey_patch['read'].add(cls)


io_registry = IORegistry()
sniff = io_registry.sniff
read = io_registry.read
write = io_registry.write

def create_format(name, encoding=None, newline=None):
    format = Format(name, encoding=encoding, newline=newline)
    io_registry.add_format(format)
    return format
