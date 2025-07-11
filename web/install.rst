Installing scikit-bio
=====================

Scikit-bio can be installed via `PyPI <https://pypi.org/>`_, `Conda <https://docs.conda.io/>`_, system package managers, from source, or from pre-compiled wheels. The PyPI and Conda methods are the most commonly used.


System requirements
-------------------

Platform
^^^^^^^^

Scikit-bio runs on Linux, macOS, and Windows. It is currently available for the x86_64 and ARM64 architectures, although ARM64 support is not yet available for the Windows operating system.

Python
^^^^^^

scikit-bio requires `Python <https://www.python.org/>`_ 3.9 or later installed in your system. See the `Python version support`_ section for more details.


Environment based
-----------------
Conda
^^^^^

The recommended way to install scikit-bio is via the Conda package manager. The latest release of scikit-bio is distributed via the `conda-forge <https://conda-forge.org/>`_ channel. You can install it via the following command::

    conda install -c conda-forge scikit-bio

Other channels such as anaconda and bioconda also host scikit-bio, which however may or may not be the up-to-date version.


PyPI
^^^^

Alternatively, the latest release of scikit-bio can be installed from PyPI::

    pip install scikit-bio


System package managers
-----------------------

Scikit-bio is available as third-party packages from software repositories for multiple Linux/BSD distributions. However, these packages may or may not be the latest version. **The scikit-bio development team is not involved in the maintenance of these packages**.

For example, users of Debian-based Linux distributions (such as Ubuntu and Linux Mint) may install scikit-bio using::

    sudo apt install python3-skbio python-skbio-doc

Users of Arch Linux or variants (such as Manjaro) may install scikit-bio from AUR::

    yay -S python-scikit-bio


Pre-compiled wheels
-------------------

Starting with version 0.7.0, scikit-bio now provides `pre-compiled wheels <https://pypi.org/project/scikit-bio/#files>`_ for each release. To install from a wheel file, download the appropriate file for your platform and run::

    pip install your_wheel_file.whl


Development version
-------------------

Scikit-bio is undergoing expansion, with many new features being introduced. You are welcome to try these features by installing the current development version from our `GitHub repo <https://github.com/scikit-bio/scikit-bio>`_::

    pip install git+https://github.com/scikit-bio/scikit-bio.git

Alternatively, you may download the repository, extract, and execute::

    pip install .

However, be cautious that the new functionality may not be stable and could be changed in the next formal release. It is not recommended to deploy the development version in a production environment.


Verifying the installation
--------------------------

After installing scikit-bio, verify the installation by running the following in a Python shell or script::

    import skbio
    print(skbio.__version__)

This should print the installed version of scikit-bio without errors.

For a more robust verification of the installation you may run the scikit-bio's unit tests in the environment where scikit-bio is installed. This will require the `pytest <https://github.com/pytest-dev/pytest>`_ package to be installed::

    python -m skbio.test


Acceleration
------------

The new `scikit-bio-binaries <https://github.com/scikit-bio/scikit-bio-binaries>`_ package is a separate package from scikit-bio, written in C++, which when installed in the same environment as scikit-bio will dramatically increase performance of select functions. Installation of scikit-bio-binaries is currently available with conda::

    conda install -c conda-forge scikit-bio-binaries


Python version support
----------------------

Scikit-bio's policy is to support Python versions until they reach official end-of-life, which typically provides 4-5 years of support per version. This approach accounts for the prevalence of legacy and under-maintained software packages in bioinformatics workflows, where users often cannot update their environments frequently. The goal is to strike an appropriate balance between encouraging modernization and maintaining compatibility for our user base, while still providing predictable deprecation timelines that researchers can plan around.
