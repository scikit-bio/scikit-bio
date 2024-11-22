# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from functools import wraps

from ._exception import OverrideError
from ._warning import _warn_renamed


def aliased(alias, deprecated=False):
    """Specify an alias for a function.

    Parameters
    ----------
    alias : str
        Alias name of the function.
    deprecated : bool or str, optional
        Whether to display a deprecation warning when the alias is called. Can also
        specify a version number indicating when the alias will be removed.

    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            _warn_renamed(func, alias, deprecated)
            return func(*args, **kwargs)

        # get parent module
        module = func.__module__
        mpath, _, mname = module.rpartition(".")
        if mname.startswith("_"):
            module = mpath
        parent = sys.modules[module]

        # get parent class (if applicable)
        # https://stackoverflow.com/questions/3589311/
        # qualname = func.__qualname__.split(".<locals>", 1)[0]
        # for level in qualname.split(".")[:-1]:
        #     parent = getattr(parent, level)

        target = wrapper if deprecated else func
        setattr(parent, alias, target)

        # add alias to docstring
        msg = f"    Alias: ``{alias}``{' (deprecated)' if deprecated else ''}.\n"
        if (doc := func.__doc__) is None:
            func.__doc__ = msg
        elif pos := doc.find("\n\n") + 1:
            func.__doc__ = doc[:pos] + "\n" + msg + doc[pos:]
        else:
            func.__doc__ += "\n" + msg

        return func

    return decorator


def meth_alias(alias_name, deprecated=False):
    def decorator(func):
        func._alias = alias_name

        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return func

    return decorator


def register_aliases(cls):
    toadd = []
    for meth in cls.__dict__.values():
        if hasattr(meth, "_alias"):
            toadd.append(meth)
    for meth in toadd:
        alias = meth._alias

        @wraps(meth)
        def wrapper(*args, **kwargs):
            _warn_renamed(meth, alias, "0.9.5")
            return meth(*args, **kwargs)

        setattr(cls, alias, wrapper)
        delattr(meth, "_alias")
    return cls


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
