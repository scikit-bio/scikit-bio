# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

TARGET ?= skbio

ifeq ($(WITH_COVERAGE), TRUE)
	COVERAGE = coverage=True
else
	COVERAGE = coverage=False
endif

test:
	python -c "import skbio; skbio.test($(COVERAGE))"
	pep8 skbio setup.py checklist.py
	flake8 skbio setup.py checklist.py
	./checklist.py
