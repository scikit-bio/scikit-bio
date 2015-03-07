# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(VERBOSE), TRUE)
	TEST_OPTIONS = verbose=True
endif

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = coverage run -m skbio.__init__
else
	TEST_COMMAND = python -c "import skbio; skbio.test($(TEST_OPTIONS))"
endif

test:
	$(TEST_COMMAND)
	pep8 skbio setup.py checklist.py
	flake8 skbio setup.py checklist.py
	./checklist.py
