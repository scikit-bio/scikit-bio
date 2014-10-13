from __future__ import absolute_import, division, print_function
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
import types
import traceback

from future.builtins import zip

from . import (UnrecognizedFormatError, InvalidRegistrationError,
               DuplicateRegistrationError, ArgumentOverrideWarning,
               FormatIdentificationWarning)
from .util import open_file, open_files, _is_string_or_bytes
from skbio.util import flatten

_formats = {}
_sniffers = {}
_aliases = {}
_empty_file_format = '<emptyfile>'


def _compound_format(fmts):
    return ', '.join(fmts)


def _factor_format(fmt):
    if _is_string_or_bytes(fmt):
        return [f.strip() for f in fmt.split(',')]
    return fmt


def _format_len(fmt):
    return len(_factor_format(fmt))


def _normalize_format(fmt):
    """Return normalized format string, is_compound format."""
    if _is_string_or_bytes(fmt):
        return _compound_format(sorted(
                                _factor_format(fmt.lower()))), ',' in fmt
    else:
        return _compound_format(sorted([_normalize_format(f)[0] for f in
                                fmt])), True


def _is_iter_list(x):
    return (hasattr(x, '__iter__') and not hasattr(x, 'read') and
            not _is_string_or_bytes(x))


def _setup_kwargs(kws):
    if len(kws) == 1:
        return kws[0]
    kwargs = {}
    for key in flatten([k.keys() for k in kws]):
        kwarg = []
        for kw in kws:
            kwarg.append(kw.get(key, None))
        kwargs[key] = kwarg
    return kwargs


def _override_kwarg(kw, key, value, warn_user):
    if key in kw and warn_user and kw[key] != value:
        warn('Best guess was: %s=%s, continuing with user supplied: %s' % (
            key, str(kw[key]), str(value)
        ), ArgumentOverrideWarning)
    kw[key] = value


def _override_kwargs(kw, fmt_kw, fmt_len, warn_user):
    for key in kw:
        if fmt_len > 1 and (not _is_iter_list(kw[key]) or
                            len(kw[key]) != fmt_len):
            _override_kwarg(fmt_kw, key, [kw[key]] * fmt_len, warn_user)
        else:
            _override_kwarg(fmt_kw, key, kw[key], warn_user)
    return fmt_kw


def register_sniffer(format):
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
    fmt, is_compound = _normalize_format(format)
    if is_compound:
        raise ValueError("'register_sniffer' cannot be applied to compound "
                         "formats.")

    def decorator(sniffer):
        if fmt in _sniffers:
            raise DuplicateRegistrationError(msg="'%s' already has a sniffer."
                                             % fmt)

        def wrapped_sniffer(fp, mode='U', **kwargs):
            with open_file(fp, mode) as fh:
                #orig_pos = fh.tell()
                #fh.seek(0)
                try:
                    return sniffer(fh, **kwargs)
                except Exception:
                    warn("'%s' has encountered a problem.\n"
                         "Please send the following to our issue tracker at\n"
                         "https://github.com/biocore/scikit-bio/issues\n\n"
                         "%s" % (sniffer.__name__, traceback.format_exc()),
                         FormatIdentificationWarning)
                    return False, {}
                finally:
                    #fh.seek(orig_pos)
                    pass

        wrapped_sniffer.__doc__ = sniffer.__doc__
        wrapped_sniffer.__name__ = sniffer.__name__

        _sniffers[fmt] = wrapped_sniffer
        return wrapped_sniffer
    return decorator


def register_reader(format, cls=None):
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

    In the event of a compound format (`['format1', 'format2']`) filehandles
    will be unrolled in the same order as the format and ALL kwarg arguments
    will be passed as tuples in the same order as the format. i.e.
    ``def format1_format2_to_generator(fmt1_fh, fmt2_fh, some_arg=(1, 2)):``

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
    fmt, is_compound = _normalize_format(format)

    def decorator(reader):
        format_class = _formats.setdefault(fmt, {}).setdefault(cls, {})

        if 'reader' in format_class:
            raise DuplicateRegistrationError('reader', fmt, cls)

        # We wrap the reader so that basic file handling can be managed
        # externally from the business logic.
        if cls is None:
            # In the case of generators, we need to make sure that the
            # filehandle is not mutated by calling `next` on the generator.
            # This allows for the odd edge case of creating several generators
            # from the same filehandle without having them inadvertently
            # interact via the filehandle's position.
            def wrapped_reader(fp, mode='U', mutate_fh=False, **kwargs):
                if not _is_iter_list(fp):
                    fp = [fp]

                with open_files(fp, mode) as fhs:
                    try:
                        generator = reader(*fhs, **kwargs)
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

                    # If the user has permitted us to mutate the file, or if
                    # we know that the registry 'owns' the filehandle then we
                    # don't need to worry about mutations. The caveat is if it
                    # is a compound format, in which case we are no longer
                    # certain who owns which file.
                    if mutate_fh or (not is_compound and
                                     _is_string_or_bytes(fp[0])):
                        while True:
                            yield next(generator)

                    else:
                        # Preserve the original filehandle positions at each
                        # call for `next`.
                        #orig_positions = [fh.tell() for fh in fhs]
                        #read_positions = orig_positions
                        try:
                            while True:
                                #orig_positions = [fh.tell() for fh in fhs]

                                #for fh, pos in zip(fhs, read_positions):
                                    #fh.seek(pos)
                                next_result = next(generator)
                                #read_positions = [fh.tell() for fh in fhs]

                                #for fh, pos in zip(fhs, orig_positions):
                                    #fh.seek(pos)

                                yield next_result
                        finally:
                            pass
                            #for fh, pos in zip(fhs, orig_positions):
                                #fh.seek(pos)

        else:
            # When an object is instantiated we don't need to worry about the
            # original position at every step, only at the end.
            def wrapped_reader(fp, mode='U', mutate_fh=False, **kwargs):
                if not _is_iter_list(fp):
                    fp = [fp]

                with open_files(fp, mode) as fhs:
                    try:
                        #orig_positions = [fh.tell() for fh in fhs]
                        return reader(*fhs, **kwargs)
                    finally:
                        if not mutate_fh:
                            pass
                            #for fh, pos in zip(fhs, orig_positions):
                                #fh.seek(pos)

        wrapped_reader.__doc__ = reader.__doc__
        wrapped_reader.__name__ = reader.__name__

        format_class['reader'] = wrapped_reader
        format_class['reader_args'] = _factor_format(format)
        return wrapped_reader
    return decorator


def register_writer(format, cls=None):
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

    In the event of a compound format (`['format1', 'format2']`) filehandles
    will be unrolled in the same order as the format and ALL kwarg arguments
    will be passed as tuples in the same order as the format. i.e.
    ``def gen_to_format1_format2(gen, fmt1_fh, fmt2_fh, some_arg=(1, 2)):``

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
    fmt, is_compound = _normalize_format(format)

    def decorator(writer):
        format_class = _formats.setdefault(fmt, {}).setdefault(cls, {})

        if 'writer' in format_class:
            raise DuplicateRegistrationError('writer', fmt, cls)

        # We wrap the writer so that basic file handling can be managed
        # externally from the business logic.
        def wrapped_writer(obj, fp, mode='w', **kwargs):
            if not _is_iter_list(fp):
                fp = [fp]

            with open_files(fp, mode) as fhs:
                writer(obj, *fhs, **kwargs)

        wrapped_writer.__doc__ = writer.__doc__
        wrapped_writer.__name__ = writer.__name__

        format_class['writer'] = wrapped_writer
        format_class['writer_args'] = _factor_format(format)
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
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                f = _formats[fmt][cls][name+'_args']
                formats.append(_compound_format(f))
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
    fmt, is_compound = _normalize_format(format)
    if not is_compound:
        if fmt in _sniffers:
            return _sniffers[fmt]
        return None
    else:
        # In the event of a compound format, we can generate a compound sniffer
        # by simply composing the results of each sub-sniffer.
        sniffers = []
        for f in _factor_format(format):
            sniffer = get_sniffer(f)
            if sniffer is None:
                return None
            sniffers.append(sniffer)

        def sniffer(fp, mode='U'):
            kwargs = []
            if not _is_iter_list(fp):
                raise ValueError('Must supply a list of files.')
            if len(fp) != len(sniffers):
                raise ValueError('List length (%d) must be %d.'
                                 % (len(fp), len(sniffers)))
            for f, sniffer in zip(fp, sniffers):
                is_format, fmt_kwargs = sniffer(f, mode=mode)
                if not is_format:
                    return False, {}
                kwargs.append(fmt_kwargs)
            return True, _setup_kwargs(kwargs)

        return sniffer


def get_reader(format, cls=None):
    """Return a reader for a format.

    Parameters
    ----------
    format : str or iterable of str
        A registered format string or compound format.
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
    fmt, is_compound = _normalize_format(format)
    composition = _factor_format(format)

    reader, original_format_order = _rw_getter('reader', fmt, cls)
    if reader is None:
        return None

    if not is_compound or original_format_order == composition:
        return reader

    # Time to generate a flip on the fly! :rimshot:
    def generated_reader(fp, **kwargs):
        if len(fp) != len(original_format_order):
            raise ValueError('List length (%d) must be %d.'
                             % (len(fp), len(original_format_order)))
        mapped_fp = [None for f in fp]
        for i, f in enumerate(original_format_order):
            mapped_fp[i] = fp[composition.index(f)]
        return reader(mapped_fp, **kwargs)

    generated_reader.__name__ = 'flip_of_' + reader.__name__
    # Generated readers do not get a copy of the docstring because
    # their argument order is not the same.

    return generated_reader


def get_writer(format, cls=None):
    """Return a writer for a format.

    Parameters
    ----------
    format : str or iterable of str
        A registered format string or compound format.
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
    fmt, is_compound = _normalize_format(format)
    composition = _factor_format(format)

    writer, original_format_order = _rw_getter('writer', fmt, cls)
    if writer is None:
        return None

    if not is_compound or original_format_order == composition:
        return writer

    def generated_writer(obj, fp, **kwargs):
        if len(fp) != len(original_format_order):
            raise ValueError('List length (%d) must be %d.'
                             % (len(fp), len(original_format_order)))
        mapped_fp = [None for f in fp]
        for i, f in enumerate(original_format_order):
            mapped_fp[i] = fp[composition.index(f)]
        return writer(obj, mapped_fp, **kwargs)

    generated_writer.__name__ = 'flip_of_' + writer.__name__
    # Generated writers do not get a copy of the docstring because
    # their argument order is not the same.

    return generated_writer


def _rw_getter(name, fmt, cls):
    if fmt in _formats:
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                return (_formats[fmt][cls][name],
                        _formats[fmt][cls][name+"_args"])
    return None, None


def sniff(fp, cls=None, mode='U'):
    """Attempt to guess the format of a file and return format str and kwargs.

    Parameters
    ----------
    fp : filepath or filehandle
        The provided file to guess the format of. Filepaths are automatically
        closed; filehandles are the responsibility of the caller.
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

    """
    if not _is_iter_list(fp):
        fp = [fp]
    factored_format = []
    kwargs = []
    for f in fp:
        possibles = []
        for fmt in _sniffers:
            if cls is not None and fmt != _empty_file_format and (
                    fmt not in _formats or cls not in _formats[fmt]):
                continue
            format_sniffer = _sniffers[fmt]
            is_format, fmt_kwargs = format_sniffer(f, mode=mode)
            if is_format:
                possibles.append(fmt)
                kwargs.append(fmt_kwargs)

        if not possibles:
            raise UnrecognizedFormatError("Cannot guess the format for %s."
                                          % str(f))
        if len(possibles) > 1:
            raise UnrecognizedFormatError("File format is ambiguous, may be"
                                          " one of %s." % str(possibles))

        factored_format.append(possibles[0])
    return _compound_format(factored_format), _setup_kwargs(kwargs)


def read(fp, format=None, into=None, verify=True, mode='U', **kwargs):
    """Read a supported skbio file format into an instance or a generator.

    This function is able to reference and execute all *registered* read
    operations in skbio.

    Parameters
    ----------
    fp : filepath, filehandle, or iterable of either
        The location to read the given `format` `into`. Filepaths are
        automatically closed when read; filehandles are the responsibility
        of the caller. In the case of a generator, a filepath will be closed
        when ``StopIteration`` is raised; filehandles are still the
        responsibility of the caller. If `format` is a compound format, then
        `fp` **must** be an iterable of the same length as the compound format.
    format : str or iterable of str, optional
        The format must be a format name with a reader for the given
        `into` class. In the case of compound formats, any order of the simple
        formats will work. If a `format` is not provided or is None, all
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

    """
    if format is None and into is None:
        raise ValueError("`format` and `into` cannot both be None.")

    if format is None:
        format, fmt_kwargs = sniff(fp, cls=into, mode=mode)
        kwargs = _override_kwargs(kwargs, fmt_kwargs, _format_len(format),
                                  verify)
    elif verify:
        sniffer = get_sniffer(format)
        if sniffer is not None:
            is_format, fmt_kwargs = sniffer(fp)
            if not is_format:
                warn("%s could not be positively identified as %s file." %
                     (str(fp), format),
                     FormatIdentificationWarning)
            else:
                kwargs = _override_kwargs(kwargs, fmt_kwargs,
                                          _format_len(format), True)

    reader = get_reader(format, into)
    if reader is None:
        raise UnrecognizedFormatError("Cannot read %s into %s, no reader "
                                      "found." % (format, into.__name__
                                                  if into is not None
                                                  else 'generator'))
    return reader(fp, mode=mode, **kwargs)


def write(obj, format, into, mode='w', **kwargs):
    """Write a supported skbio file format from an instance or a generator.

    This function is able to reference and execute all *registered* write
    operations in skbio.

    Parameters
    ----------
    obj : object
        The object must have a registered writer for a provided `format`.
    format : str or iterable of str
        The format must be a registered format name with a writer for the given
        `obj`. In the case of compound formats, any order of the simple
        formats will work.
    into : filepath, filehandle or iterable of either
        The location to write the given `format` from `obj` into. Filepaths are
        automatically closed when written; filehandles are the responsibility
        of the caller. If `format` is a compound format, then `into` **must**
        be an iterable of the same length as the compound format.
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

    writer(obj, into, mode=mode, **kwargs)


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
        fmts = _factor_format(fmt)
        lines += ("skbio.io." + fmt for fmt in fmts)
    return '\n'.join(lines)


def _formats_for_docs(formats):
    lines = []
    for fmt in formats:
        fmts = _factor_format(fmt)
        if len(fmts) > 1:
            lines.append("- ``%s`` (%s)" % (
                fmts, ', '.join(":mod:`skbio.io.%s`" % s for s in fmts)))
        else:
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
fp : filepath, filehandle, or iterable of either
    The location to read the given `format`. Filepaths are
    automatically closed when read; filehandles are the
    responsibility of the caller. If `format` is a compound format,
    then `fp` **must** be an iterable of the same length as the
    compound format.
format : str or iterable of str, optional
    The format must be a format name with a reader for ``%s``.
    In the case of compound formats, any order of the simple
    formats will work as long as `fp` follows the same order.
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
