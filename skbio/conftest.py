# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

def pytest_configure(config):
    """Add configuration options for the doctests

    Conventionally this configuration would go in the ini-options file, but we
    need to be able to run the tests in an installed package which would need
    a configuration file for every installation.
    """
    config.addinivalue_line("doctest_optionflags", "NORMALIZE_WHITESPACE")
    config.addinivalue_line("doctest_optionflags", "IGNORE_EXCEPTION_DETAIL")
