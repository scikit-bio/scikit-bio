# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = coverage run --rcfile ../.coveragerc -m skbio.test && coverage report --rcfile ../.coveragerc
else
	TEST_COMMAND = python -m skbio.test
endif

.PHONY: doc web lint test dev install

doc:
	$(MAKE) -C doc clean html

web:
	$(MAKE) -C web clean html

clean:
	$(MAKE) -C doc clean
	$(MAKE) -C web clean
	rm -rf build dist scikit_bio.egg-info

lint:
	ruff check skbio setup.py checklist.py
	./checklist.py
	check-manifest

# cd into a directory that is different from scikit-bio root directory to
# simulate a user's install and testing of scikit-bio. Running from the root
# directory will find the `skbio` subpackage (not necessarily the installed
# one!) because cwd is considered in Python's search path. It is important to
# simulate a user's install/test process this way to find package data that did
# not install correctly (for example).
test:
	cd ci && $(TEST_COMMAND)

install:
	pip install .

dev:
	pip install -e .
