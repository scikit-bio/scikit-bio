from __future__ import division, absolute_import, print_function

from .numpydoc import setup

# This is so that our docs build.
def _closure():
    from skbio.util import classproperty
    def __get__(self, cls, owner):
        return self

    classproperty.__get__ = __get__

_closure()
