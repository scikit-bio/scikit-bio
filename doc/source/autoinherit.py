# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from sphinx.util.inspect import safe_getattr


def trace_inherit(name, real_name):
    """Trace the base class from which a member is inherited.

    Parameters
    ----------
    name : str
        Public name of the member.
    real_name : str
        True import path of the member.

    Returns
    -------
    str
        Import path of the member from the base class.

    """
    # `real_name` is the full path to the member, such as:
    #   skbio.sequence.DNA.frequencies
    # It can be split into:
    #   module: skbio.sequence
    #   class: DNA
    #   member: frequencies
    # `member` should be identical to `name`, therefore the `name` parameter
    # is actually redundant in this function.
    try:
        mod_name, cls_name, _ = real_name.rsplit(".", 2)
    except ValueError:
        return

    # Import the class which has the member.
    try:
        module = __import__(mod_name, fromlist=[cls_name])
    except ImportError:
        return
    try:
        cls = safe_getattr(module, cls_name, None)
    except AttributeError:
        return

    # Make sure class has the member.
    # If the member is inherited, `hasattr` returns True but `cls.__dict__`
    # does not have it.
    if not hasattr(cls, name):
        return

    # Trace back the method resolution order (MRO) of the class to identify
    # the base class that has the same member.
    # If the member is owned by the current class (i.e., not inherited), it
    # will return the current class.
    origin = None
    for base in cls.__mro__:
        if name in base.__dict__:
            if base.__module__.startswith("skbio."):
                origin = base
            break
    if not origin:
        return

    # Get fully qualified name of the base class, such as:
    #   skbio.sequence.Sequence
    # In scikit-bio, classes are implemented in files starting with `_` (e.g.,
    # `_base.py`) but are imported in `__init__.py` of the module. Therefore,
    # the following code converts the name below into the above:
    #   skbio.sequence._sequence.Sequence
    # However, if in the future this structure changes, the code needs to be
    # modified.
    path = origin.__module__.split(".")
    while path and path[-1].startswith("_"):
        del path[-1]

    # Eventually, it returns:
    #   skbio.sequence.Sequence.frequencies
    return ".".join(path + [origin.__name__, name])
