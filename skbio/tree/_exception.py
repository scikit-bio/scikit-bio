# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class TreeError(Exception):
    """General tree error"""
    pass


class NoLengthError(TreeError):
    """Missing length when expected"""
    pass


class DuplicateNodeError(TreeError):
    """Duplicate nodes with identical names"""
    pass


class MissingNodeError(TreeError):
    """Expecting a node"""
    pass


class NoParentError(MissingNodeError):
    """Missing a parent"""
    pass
