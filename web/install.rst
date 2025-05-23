Install scikit-bio
==================


Python
------

scikit-bio requires `Python <https://www.python.org/>`_ 3.9 or later installed in your system.


Conda
-----

The recommended way to install scikit-bio is via the `Conda <https://docs.conda.io/>`_ package manager. The latest release of scikit-bio is distributed via the `conda-forge <https://conda-forge.org/>`_ channel. You can install it via the following command::

    conda install -c conda-forge scikit-bio

Other channels such as anaconda and bioconda also host scikit-bio, which however may or may not be the up-to-date version.


PyPI
----

Alternatively, the latest release of scikit-bio can be installed from `PyPI <https://pypi.org/>`_::

    pip install scikit-bio


Third-party
-----------

scikit-bio is available as third-party packages from software repositories for multiple Linux/BSD distributions. However, these packages may or may not be the latest version. The scikit-bio development team is not involved in the maintenance of these packages.

For example, users of Debian-based Linux distributions (such as Ubuntu and Linux Mint) may install scikit-bio using::

    sudo apt install python3-skbio python-skbio-doc

Users of Arch Linux or variants (such as Manjaro) may install scikit-bio from AUR::

    yay -S python-scikit-bio


Nightly build
-------------

scikit-bio is undergoing expansion, with many new features being introduced. You are welcome to try these features by installing the current development version from our `GitHub repo <https://github.com/scikit-bio/scikit-bio>`_.

    pip install git+https://github.com/scikit-bio/scikit-bio.git

Alternatively, you may download the repository, extract, and execute::

    python setup.py install

However, be cautious that the new functionality may not be stable and could be changed in the next formal release. It is not recommended to deploy the development version in a production environment.


.. Test
.. ----

.. You can verify your installation by running the scikit-bio unit tests (this requires `pytest` installed)::

..     python -m skbio.test

.. If the installation was successful and all features of scikit-bio work as intended, the test will report only passes (and warnings), but no failures.
