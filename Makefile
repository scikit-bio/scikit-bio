# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = python setup.py test -a "--cov=skbio"
else
	TEST_COMMAND = python setup.py test
endif

test:
	$(TEST_COMMAND)
	flake8 skbio setup.py checklist.py
	./checklist.py
	check-manifest
