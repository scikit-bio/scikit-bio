# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_COVERAGE), TRUE)
	# TODO: actually add coverage back in
	TEST_COMMAND = python setup.py test
else
	TEST_COMMAND = python setup.py test
endif

test:
	$(TEST_COMMAND)
	flake8 skbio setup.py checklist.py
	./checklist.py
	check-manifest
