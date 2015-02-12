# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

TARGET ?= skbio

ifeq ($(WITH_COVERAGE), TRUE)
	COVERAGE = --with-coverage
endif

ifeq ($(WITH_DOCTEST), FALSE)
	DOCTEST =
else
	DOCTEST = --with-doctest
endif

test:
	nosetests $(TARGET) $(DOCTEST) $(COVERAGE) --ignore DO_NOT_IGNORE_ANYTHING
	pep8 $(TARGET) setup.py checklist.py
	flake8 $(TARGET) setup.py checklist.py
	./checklist.py
