# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import wraps
from inspect import signature

from ._exception import OverrideError
from ._warning import _warn_renamed


# def func_alias(alias, deprecated=False):
#     """Create an alias for a function.

#     Parameters
#     ----------
#     alias : str
#         Alias name of the function.
#     deprecated : bool or str, optional
#         Whether to display a deprecation warning when the alias is called. Can also
#         specify a version number indicating when the alias will be removed.

#     """

#     def decorator(func):
#         @wraps(func)
#         def wrapper(*args, **kwargs):
#             _warn_renamed(func, alias, deprecated)
#             return func(*args, **kwargs)

#         # get parent module
#         module = func.__module__
#         mpath, _, mname = module.rpartition(".")
#         if mname.startswith("_"):
#             module = mpath
#         parent = sys.modules[module]

#         # get parent class, if applicable
#         # see: https://stackoverflow.com/questions/3589311/
#         # this code however doesn't work as it triggers circular import
#         qualname = func.__qualname__.split(".<locals>", 1)[0]
#         for level in qualname.split(".")[:-1]:
#             parent = getattr(parent, level)

#         target = wrapper if deprecated else func
#         setattr(parent, alias, target)

#         # add alias to docstring
#         msg = f"    Alias: ``{alias}``{' (deprecated)' if deprecated else ''}.\n"
#         if (doc := func.__doc__) is None:
#             func.__doc__ = msg
#         elif pos := doc.find("\n\n") + 1:
#             func.__doc__ = doc[:pos] + "\n" + msg + doc[pos:]
#         else:
#             func.__doc__ += "\n" + msg

#         return func

#     return decorator


def alias_doc(func, name, since=None, until=None, warn=False, params={}):
    """Modify the docstring of a function to indicate the alias status."""
    if not since:
        msg = f"Alias: ``{name}``"
    else:
        msg = (
            f".. versionchanged:: {since} "
            f"Renamed from ``{name}``. The old name is kept as an alias"
        )
        if warn:
            msg += " but is deprecated"
            if until:
                msg += f" and will be removed in {until}"
        msg += "."

    func.__doc__ = alias_into_doc(func.__doc__, msg)


def alias_into_doc(doc, msg):
    """Insert text into the docstring of a function."""
    # no docstring: message becomes the entire docstring.
    if doc is None:
        return msg + "\n"

    lines = doc.splitlines()
    n = len(lines)

    # find indentation size, which is important for the docstring to be rendered
    indent = 0
    for line in lines[1:]:
        if line:
            indent = len(line) - len(line.lstrip())
            break

    indent = " " * indent

    # docstring has at least two paragraphs: insert message after the first paragraph
    for i in range(1, n):
        if not lines[i]:
            return "\n".join(lines[:i] + ["", indent + msg] + lines[i:])

    # docstring has only one paragraph: append message to the end of docstring
    return doc + "\n" + indent + msg + "\n"


def aliased(name, since=None, until=None, warn=False, params={}):
    """Create an alias for a function or method.

    Parameters
    ----------
    name : str
        Alias name of the function or method.
    since : str, optional
        Version when alias was created.
    until : str, optional
        Version when alias will be removed.
    warn : bool, optional
        Raise a deprecation warning when alias is called.
    params : dict, optional
        Aliases of parameters. Each key is a parameter, and its value consists of
        name, since, until and warn.

    See Also
    --------
    add_aliases

    Notes
    -----
    This is a decorator that can be applied to a function or method to indicate its
    alias status.

    """

    def decorator(func):
        func._alias = (name, since, until, warn, params)
        return func

    return decorator


def add_aliases(cls):
    """Register aliases of members of a module or class.

    For each member that has an alias, create a wrapper named as the alias, and add a
    "skip" flag which instructs Sphinx autodoc to skip it when rendering the
    documentation.

    See Also
    --------
    aliased

    """
    toadd = []
    for func in cls.__dict__.values():
        if hasattr(func, "_alias"):
            toadd.append(func)

    for func in toadd:
        name, since, until, warn, params = func._alias

        @wraps(func)
        def wrapper(*args, **kwargs):
            if warn:
                _warn_renamed(func, name, since, until)
            return func(*args, **kwargs)

        wrapper._skip = True
        delattr(wrapper, "_alias")
        setattr(cls, name, wrapper)

        alias_doc(func, name, since, until, warn, params)

    return cls


def params_aliased(params=[]):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for key, alias, since, until, warn in params:
                if alias in kwargs:
                    kwargs[key] = kwargs.pop(alias)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def overrides(interface_class):
    """Indicate that a member is being overridden from a specific parent class.

    Decorator for class-level members. Used to indicate that a member is being
    overridden from a specific parent class. If the member does not have a docstring,
    it will pull one from the parent class. When chaining decorators, this should be
    first as it is relatively nondestructive.

    Parameters
    ----------
    interface_class : class
        The class which has a member overridden by the decorated member.

    Returns
    -------
    function
        The function is not changed or replaced.

    Raises
    ------
    OverrideError
        If the `interface_class` does not possess a member of the same name
        as the decorated member.

    """
    # Adapted from http://stackoverflow.com/a/8313042/579416

    def overrider(method):
        if method.__name__ not in dir(interface_class):
            raise OverrideError(
                f"{method.__name__} is not present in parent "
                f"class: {interface_class.__name__}."
            )
        backup = classproperty.__get__
        classproperty.__get__ = lambda x, y, z: x
        if method.__doc__ is None:
            method.__doc__ = getattr(interface_class, method.__name__).__doc__
        classproperty.__get__ = backup
        return method

    return overrider


class classproperty(property):
    """Decorator for class-level properties.

    Supports read access only. The property will be read-only within an
    instance. However, the property can always be redefined on the class, since
    Python classes are mutable.

    Parameters
    ----------
    func : function
        Method to make a class property.

    Returns
    -------
    property
        Decorated method.

    Raises
    ------
    AttributeError
        If the property is set on an instance.

    """

    def __init__(self, func):
        name = func.__name__
        doc = func.__doc__
        super(classproperty, self).__init__(classmethod(func))
        self.__name__ = name
        self.__doc__ = doc

    def __get__(self, cls, owner):
        return self.fget.__get__(None, owner)()

    def __set__(self, obj, value):
        raise AttributeError("can't set attribute")


class classonlymethod(classmethod):
    """Just like `classmethod`, but it can't be called on an instance."""

    def __get__(self, obj, cls=None):
        if obj is not None:
            raise TypeError(
                f"Class-only method called on an instance. Use "
                f"{cls.__name__}.{self.__func__.__name__} "
                "instead."
            )
        return super().__get__(obj, cls)
