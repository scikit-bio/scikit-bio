# Releasing a new version of scikit-bio

## Introduction

Releasing a piece of software can simultaneously be an invigorating, intimidating, horrifying, and cathartic experience. This guide aims to make the release process as smooth as possible, though this goal may be impossible to fully meet.

To illustrate examples of commands you might run, let's assume that the current version is 1.2.3-dev and we want to release version 1.2.4.

**Note:** the following commands assume you are in the top-level directory of the scikit-bio repository unless otherwise noted.

**Tip:** It can be efficient to have the help of a couple other devs as some steps can be run in parallel. It's also useful to have a variety of platforms/environments to test on during the release process, so find friends that are Linux and Mac users!

## Prepping the release

1. Ensure the scikit-bio build is passing on Travis.

2. Update the version strings (1.2.3-dev) to the new version (1.2.4). There should only be two places this needs to be done: ``setup.py`` and ``skbio.__init__.py``. It's a good idea to ``grep`` for the current version string just to be safe:

    grep -ir '1\.2\.3-dev' *

Update ``CHANGELOG.md`` to include descriptions of the changes that made it into this release. Be sure to update the heading to include the new version (1.2.4) and the date of the release. Use the existing structure in the file as a template/guide.

Submit a pull request with these changes and let Travis run the tests. In the meantime, build the documentation locally:

    cd doc
    make clean && make html

**Note:** you will need to **fully install** (including built extensions) the exact version of scikit-bio that you are editing so that Sphinx will pull docstrings from the correct version of the code. **Make sure the version of scikit-bio that is imported by ``import skbio`` is the correct one!**

Switch to the ``gh-pages`` branch of the repository; this is where the website and built documentation is stored.

Remove everything from ``docs/latest/``:

    git rm -rf docs/latest/*

Create a directory for the new version of the docs and recreate the ``latest/`` directory:

    mkdir docs/1.2.4
    mkdir docs/latest

Copy over the built documentation to both ``docs/1.2.4/`` and ``docs/latest``:

    cp -r <path to skbio repo>/doc/build/html/* docs/1.2.4/
    cp -r <path to skbio repo>/doc/build/html/* docs/latest/

Add a new list item to ``index.html`` to link to ``docs/1.2.4/index.html``.

Commit and push (either directly or as a pull request) to have the website updated.

If the tests passed on Travis, merge the pull request to update the version strings.

From the [scikit-bio GitHub page](https://github.com/biocore/scikit-bio), click on the releases tab and draft a new release. Use the version number for the tag name (1.2.4) and have it create the tag against the latest master branch. Fill in a release title that is consistent with the other release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea).

Once the release is created on GitHub, it's a good idea to test out the release tarball before publishing to PyPI.

Create a new virtualenv. Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory. Try installing the release and run the tests:

    pip install numpy
    pip install .
    cd
    nosetests --with-doctest skbio

During this process (it can take awhile to install all of scikit-bio's dependencies), submit a pull request to update the version strings from 1.2.4 to 1.2.4-dev. Use the same strategy described above to update the version strings. Update ``CHANGELOG.md`` to include a new section for 1.2.4-dev (there won't be any changes to note here yet). **Do not merge this pull request yet.**

Assuming the GitHub release tarball correctly installed and passed its tests, we're now ready to simulate creating the source distribution (``sdist``) that will be published to PyPI.

**Important:** Check ``MANIFEST.in`` to ensure that the files and directories it references still exist. Some may have been removed, renamed, or there may be new files/dirs that need to be included in the release. This step in the release process has caused most of the headaches; don't neglect ``MANIFEST.in``!

Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory. Build a source distribution:

    python setup.py sdist

Create a new virtualenv and run:

    cd
    pip install numpy
    pip install <path to extracted scikit-bio release>/dist/scikit-bio-0.1.3.tar.gz
    nosetests --with-doctest skbio

If everything goes well, it is finally time to push the release to PyPI:

    python setup.py sdist upload

You must have the proper login credentials to add a release to PyPI.

Once the release is available on PyPI, do a final round of testing. Create a new virtualenv and run:

    cd
    pip install numpy
    pip install scikit-bio
    nosetests --with-doctest skbio

If this succeeds, the release appears to be a success!

Merge the latest pull request to update version strings to 1.2.4-dev. Close the release milestone on the GitHub issue tracker. Send an email to the skbio developers and users groups (and anyone else who might be interested), tweet about the release, and rejoice!
