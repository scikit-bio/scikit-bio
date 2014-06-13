# Releasing a new version of scikit-bio

## Introduction

Releasing a piece of software can simultaneously be an invigorating, intimidating, horrifying, and cathartic experience. This guide aims to make the release process as smooth as possible.

To illustrate examples of commands you might run, let's assume that the current version is 1.2.3-dev and we want to release version 1.2.4.

**Note:** the following commands assume you are in the top-level directory of the scikit-bio repository unless otherwise noted.

**Tip:** It can be efficient to have the help of a couple other devs, as some steps can be run in parallel. It's also useful to have a variety of platforms/environments to test on during the release process, so find friends that are Linux and Mac users!

## Prep the release

1. Ensure the Travis build is passing against the latest master branch.

2. Update the version strings (1.2.3-dev) to the new version (1.2.4). There should only be two places this needs to be done: ``setup.py`` and ``skbio/__init__.py``. It's a good idea to ``grep`` for the current version string just to be safe:

        $ grep -ir '1\.2\.3-dev' *

3. Update ``CHANGELOG.md`` to include descriptions of the changes that made it into this release. Be sure to update the heading to include the new version (1.2.4) and the date of the release. Use the existing structure in the file as a template/guide.

4. Submit a pull request with these changes and let Travis run.

## Build the documentation

In the meantime, we can build the documentation and update the website.

**Note:** you will need to **fully install** (including built extensions) the exact version of scikit-bio that you are editing so that Sphinx will pull docstrings from the correct version of the code. **Make sure the version of scikit-bio that is imported by ``import skbio`` is the correct one!**

1. Build the documentation locally:

        cd doc
        make clean && make html

2. Switch to the ``gh-pages`` branch of the repository.

3. Remove everything from ``docs/latest/``:

        git rm -rf docs/latest/*

4. Create a directory for the new version of the docs and recreate the ``latest/`` directory:

        mkdir docs/1.2.4
        mkdir docs/latest

5. Copy over the built documentation to both ``docs/1.2.4/`` and ``docs/latest``:

        cp -r <path to skbio repo>/doc/build/html/* docs/1.2.4/
        cp -r <path to skbio repo>/doc/build/html/* docs/latest/

6. Add a new list item to ``index.html`` to link to ``docs/1.2.4/index.html``.

7. Commit and push (either directly or as a pull request) to have the website updated. This updates the live website so be sure to poke through the site to make sure things aren't broken and that version strings are correct.

8. If the tests passed on Travis, merge the pull request to update the version strings to 1.2.4.

## Tag the release

From the [scikit-bio GitHub page](https://github.com/biocore/scikit-bio), click on the releases tab and draft a new release. Use the version number for the tag name (1.2.4) and create the tag against the latest master branch. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea).

Once the release is created on GitHub, it's a good idea to test out the release tarball before publishing to PyPI.

1. Create a new virtualenv.

2. Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory.

3. Innstall the release and run the tests:

        pip install numpy
        pip install .
        cd
        nosetests --with-doctest skbio

4. During this process (it can take awhile to install all of scikit-bio's dependencies), submit a pull request to update the version strings from 1.2.4 to 1.2.4-dev. Use the same strategy described above to update the version strings. Update ``CHANGELOG.md`` to include a new section for 1.2.4-dev (there won't be any changes to note here yet). **Do not merge this pull request yet.**

## Test the source distribution

Assuming the GitHub release tarball correctly installs and passes its tests, we're now ready to test the creation of the source distribution (``sdist``) that will be published to PyPI. It is important to test the source distribution because it is created in an entirely different way than the release tarball on GitHub. Thus, there is the danger of having two different release tarballs: the one created on GitHub and the one uploaded to PyPI.

**Important:** Check ``MANIFEST.in`` to ensure that the files and directories it references still exist. Some may have been removed, renamed, or there may be new files/dirs that need to be included in the ``sdist`` release. This step in the release process has caused the most hangups; don't neglect ``MANIFEST.in``!

1. Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory.

2. Build a source distribution:

        python setup.py sdist

3. Create a new virtualenv and run:

        cd
        pip install numpy
        pip install <path to extracted scikit-bio release>/dist/scikit-bio-1.2.4.tar.gz
        nosetests --with-doctest skbio

4. If everything goes well, it is finally time to push the release to PyPI:

        python setup.py sdist upload

    You must have the proper login credentials to add a release to PyPI.

5. Once the release is available on PyPI, do a final round of testing. Create a new virtualenv and run:

        cd
        pip install numpy
        pip install scikit-bio
        nosetests --with-doctest skbio

If this succeeds, the release appears to be a success!

## Post-release cleanup

1. Merge the latest pull request to update version strings to 1.2.4-dev.

2. Close the release milestone on the GitHub issue tracker.

3. Send an email to the skbio users and developers lists, and anyone else who might be interested, e.g., lab mailing lists.

4. Tweet about the release.

5. Rejoice! no cog yay
