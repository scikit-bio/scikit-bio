from .util import open_file, _get_filehandle
from skbio.io import (FormatIdentificationError, FileFormatError,
                      DuplicateRegistrationError)
_formats = {}
_identifiers = {}


def register_identifier(fmt):
    """
    """
    def decorator(identifier):
        if fmt in _identifiers:
            raise DuplicateRegistrationError("'%s' already has an identifier."
                                             % fmt)
        _identifiers[fmt] = identifier
        return identifier
    return decorator


def register_reader(fmt, *args):
    """
    """
    return _rw_decorator('reader', fmt, *args)


def register_writer(fmt, *args):
    """
    """
    return _rw_decorator('writer', fmt, *args)


def _rw_decorator(name, fmt, *args):
    cls = None
    arg_len = len(args)
    if arg_len > 1:
        raise TypeError("register_%s takes 1 or 2 arguments (%d given)"
                        % (name, arg_len))
    if arg_len == 1:
        cls = args[0]

    def decorator(func):
        if fmt not in _formats:
            _formats[fmt] = {}
        format_dict = _formats[fmt]
        if cls not in format_dict:
            format_dict[cls] = {}
        format_class = format_dict[cls]
        if name not in format_class:
            format_class[name] = func
        else:
            raise DuplicateRegistrationError("'%s' already has a %s for %s."
                                             % (fmt, name, cls.__name__))

        return func
    return decorator


def list_read_formats(cls):
    """
    """
    return _rw_list_formats('reader', cls)


def list_write_formats(cls):
    """
    """
    return _rw_list_formats('writer', cls)


def _rw_list_formats(name, cls):
    formats = []
    for fmt in _formats:
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                formats.append(fmt)
    return formats


def get_identifier(fmt):
    """
    """
    if fmt in _identifiers:
        return _identifiers[fmt]
    return None


def get_reader(fmt, *args):
    """
    """
    return _rw_getter('reader', fmt, *args)


def get_writer(fmt, *args):
    """
    """
    return _rw_getter('writer', fmt, *args)


def _rw_getter(name, fmt, *args):
    cls = None
    arg_len = len(args)
    if arg_len > 1:
        raise TypeError("get_%s takes 1 or 2 arguments (%d given)"
                        % (name, arg_len))
    if arg_len == 1:
        cls = args[0]

    if fmt in _formats:
        if cls in _formats[fmt]:
            if name in _formats[fmt][cls]:
                return _formats[fmt][cls][name]
    return None


def guess_format(fh, cls=None):
    """
    """
    with open_file(fh, 'U') as fh:
        possibles = []
        fh.seek(0)
        for fmt in _identifiers:
            if cls is not None and (fmt not in _formats or
                                    cls not in _formats[fmt]):
                continue
            test = _identifiers[fmt]
            if test(fh):
                possibles.append(fmt)
            fh.seek(0)
        if not possibles:
            raise FormatIdentificationError("Cannot guess the format for %s."
                                            % str(fh))
        if len(possibles) > 1:
            raise FormatIdentificationError("File format is ambiguous, may be"
                                            " one of %s." % str(possibles))

        return possibles[0]


def read(fp, format=None, into=None, mode='U', *args, **kwargs):
    """
    """
    if format is None and into is None:
        raise ValueError("`format` and `into` cannot both be None.")

    fh, is_own = _get_filehandle(fp, mode)

    if format is None:
        format = guess_format(fh, cls=into)

    reader = get_reader(format, into)
    if reader is None:
        raise FileFormatError("Cannot read %s into %s, no reader found."
                              % (format, into.__name__
                                 if into is not None
                                 else 'generator'))

    if into is None:
        def wrapper_generator():
            original = reader(fh, *args, **kwargs)
            try:
                while(True):
                    yield original.next()
            finally:
                if is_own:
                    fh.close()

        return wrapper_generator()

    else:
        result = reader(fh, *args, **kwargs)
        if is_own:
            fh.close()

        return result


def write(obj, format=None, into=None, mode='w', *args, **kwargs):
    """
    """
    if format is None:
        raise ValueError("Must specify a `format` to write out as.")
    if into is None:
        raise ValueError("Must provide a filepath or filehandle for `into`")

    with open_file(into, mode) as fh:
        writer = get_writer(format, obj.__class__)
        if writer is None:
            raise FileFormatError("Cannot write %s into %s, no writer found."
                                  % (format, str(fh)))

        writer(obj, fh, *args, **kwargs)
