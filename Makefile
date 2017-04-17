# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# CONDA_ROOT is the path to the conda install directory: ${HOME}/miniconda3
# .bash_profile is reloaded in case conda was recently installed and the user
# has not updated their environment.
CONDA_ROOT := $(shell source ${HOME}/.bash_profile && conda env list 2>/dev/null | awk '$$1 == "root" {print $$NF}')

ifeq ($(WITH_COVERAGE), TRUE)
	TEST_COMMAND = COVERAGE_FILE=../.coverage coverage run \
	--rcfile ../.coveragerc -m skbio.test
else
	TEST_COMMAND = python -m skbio.test
endif

# cd into a directory that is different from scikit-bio root directory to
# simulate a user's install and testing of scikit-bio. Running from the root
# directory will find the `skbio` subpackage (not necessarily the installed
# one!) because cwd is considered in Python's search path. It is important to
# simulate a user's install/test process this way to find package data that did
# not install correctly (for example).
test:
	cd ci && $(TEST_COMMAND)
	flake8 skbio setup.py checklist.py
	./checklist.py
	check-manifest

# Automatically setup a development environment.
# If ${CONDA_ROOT}/envs/skbio exists and the 'conda' command is not found,
# prompt the user to update their environment.
# Otherwise, prompt the user to activate their environment.
# If both are active, recommend resources for getting started.
.PHONY: develop
develop: ${CONDA_ROOT}/envs/skbio
	@command -v conda >/dev/null && conda env list | \
		awk '$$1 != "skbio" && $$2 == "*" \
			{print "Run the command \x27source activate skbio\x27 to activate your development environment"; next} \
		$$1 == "skbio" && $$2 == "*" \
			{print "Looks like your environment is already setup. For more information, checkout some of these resources:\n"\
			"https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md\n"\
			"https://github.com/biocore/scikit-bio/issues\n"\
			"https://gitter.im/biocore/scikit-bio"}' || \
		echo "Run the command 'source ${HOME}/.bash_profile' or open a new terminal to finish installing 'conda'" \
		"Afterwards, run the command 'source activate skbio' to activate your development environment"; \


# See https://github.com/biocore/scikit-bio/blob/master/CONTRIBUTING.md#setting-up-a-development-environment
${CONDA_ROOT}/envs/skbio: | miniconda
	source ${HOME}/.bash_profile && \
	conda create -n skbio python=3.4 pip && \
	source activate skbio && \
	conda install --file ci/conda_requirements.txt && \
	pip install --no-deps -e .
	-$(MAKE) test

# If conda is not available, download the miniconda install script using either
# 'curl' or 'wget'. Run the install script.
.PHONY: miniconda
miniconda:
	@source ${HOME}/.bash_profile && \
	command -v conda >/dev/null || \
	{ \
		if [[ $$(uname) == "Darwin" ]]; then \
			MINICONDA_SCRIPT="Miniconda3-latest-MacOSX-x86_64.sh"; \
		else \
			MINICONDA_SCRIPT="Miniconda3-latest-Linux-x86_64.sh"; \
		fi; \
		MINICONDA_URL="https://repo.continuum.io/miniconda/$${MINICONDA_SCRIPT}"; \
		echo "The recommended development environment for contributing to scikit-bio is using Anaconda by Continuum Analytics"; \
		echo "Downloading $$MINICONDA_URL"; \
		{ command -v curl >/dev/null && curl -O -L $${MINICONDA_URL}; } || \
		{ command -v wget >/dev/null && wget $${MINICONDA_URL} -O $${MINICONDA_SCRIPT}; } && \
		bash $${MINICONDA_SCRIPT} && rm $${MINICONDA_SCRIPT}; \
	}
