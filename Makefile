# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = PYTHONSAFEPATH=1 uv run --group test coverage run --rcfile .coveragerc -m skbio.test && uv run --group test coverage report --rcfile .coveragerc
else
	TEST_COMMAND = uv run --group test python -P -m skbio.test
endif

.PHONY: doc web lint test dev install cython

doc:
	uv run --group doc $(MAKE) -C doc clean html

web:
	uv run --group doc $(MAKE) -C web clean html

clean:
	uv run --group doc $(MAKE) -C doc clean
	uv run --group doc $(MAKE) -C web clean
	rm -rf build dist scikit_bio.egg-info

lint:
	# uv run --group lint ruff check skbio setup.py checklist.py
	uv run --group lint ./checklist.py
	# uv run --group lint check-manifest

# Python 3.11+ -P / PYTHONSAFEPATH keeps the current directory off sys.path.
test:
	$(TEST_COMMAND)

install:
	uv sync

dev:
	uv sync --group dev
