# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

TARGET ?= skbio

NOSE_OPTIONS =--ignore DO_NOT_IGNORE_ANYTHING

ifeq ($(VERBOSE), TRUE)
	NOSE_OPTIONS +=-v
endif

ifeq ($(WITH_COVERAGE), TRUE)
	NOSE_OPTIONS +=--with-coverage
endif

ifneq ($(WITH_DOCTEST), FALSE)
	NOSE_OPTIONS +=--with-doctest
endif

test:
	nosetests $(NOSE_OPTIONS) $(TARGET)
	pep8 $(TARGET) setup.py checklist.py
	flake8 $(TARGET) setup.py checklist.py
	./checklist.py
