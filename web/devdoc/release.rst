Releasing a new version
=======================


Introduction
------------

This guide explains how to release a new version of scikit-bio. To illustrate examples of commands you might run, let's assume that the current version is **x.y.z-dev** and we want to release version **x.y.z**.

.. note:: The following commands assume you are in the top-level directory of the scikit-bio repository unless otherwise noted. They also assume that you have [Conda](https://conda.io) installed.

Large portions of the workflow are now automated through GitHub workflows. In order to check that things will work correctly before running the release workflow in ``release.yml``, please follow these steps.

The general workflow is:

1. Update the version number in the repository from **x.y.z-dev** to **x.y.z**.

2. Ensure documentation, wheels, and source distribution build correctly before publishing release.

3. Publish release.

4. Update version number in the repository from **x.y.z** to **a.b.c-dev**.


Prep the release
----------------

The purpose of these steps is to update the version within the package. At this point we are only updating the version number in the repository, and doing some basic checks.

1. Ensure the GitHub Actions CI build is passing against main.

2. Update the version strings (x.y.y-dev) to the new version (x.y.z). This will include ``__version__`` defined in ``skbio/__init__.py``. ``grep`` for the current version string to find all occurrences::

    grep -r 'x\.y\.y-dev' .

3. Remove any deprecated functionality that was scheduled for removal on or before this release. To find all deprecated functionality, search for `@deprecated` decorators:

    grep -r '@deprecated' .

.. note:: Any deprecated functionality that was removed in the ``### Miscellaneous`` section of ``CHANGELOG.md``.

4. Update ``CHANGELOG.md`` to include descriptions of all **user-facing** changes that made it into this release. Be sure to update the heading to include the new version (x.y.z) and the date of the release. Use the existing structure in the file as a template/guide.

5. Submit a pull request and merge when CI tests are passing.


Build documentation locally
---------------------------

This is a **pre-release check** to be performed by the developer to ensure that things run smoothly once we get to the actual release publication steps. It is a **critical step** to avoid failures once the automated workflows run. 

In an environment with **the same version of Python** specified in the ``Publish documentation`` job within the ``release.yml`` workflow, build the documentation locally using ``make doc`` and ensure that it runs without error and that everything renders correctly. 

We now use Sphinx and automated workflows to handle the documentation building and publishing process. Check Sphinx's most recent Python version support, and pin to this version in ``release.yml`` if it isn't the latest.


Test wheels and sdist fully
---------------------------

This is a **pre-release check** to be performed by the developer to ensure that things run smoothly once we get to the actual release publication steps. It is a **critical step** to avoid failures once the automated workflows run.

The ``wheels.yml`` file may be triggered manually through GitHub Actions ``workflow_dispatch`` (https://github.com/scikit-bio/scikit-bio/actions/workflows/wheels.yml). Click the ``Run workflow`` button and run it against the ``main`` branch. 

When this workflow is manually triggered, it will build wheels for all supported Python versions and all supported platforms, and run the full suite of unit tests on every wheel that it builds.
This process typically takes some time, but it is well worth the time to ensure that the wheel building process will be successful when it's time to actually publish the release.

Manually check that the correct wheels are built by downloading the artifacts from the workflow run or visually checking the output from ``cibuildwheel`` on the GitHub Actions workflow page.
If downloading, inspect the zip file containing the wheels provided by workflow for each platform. The zip file should contain all the wheels you are trying to build.

If you find that the desired wheels are not present in the zip files, or you are coming accross unexpected errors in the wheel building process, there are a few places to check:

1. Check that you have updated the ``cibuildwheel`` portion of the ``pyproject.toml`` configuration file to include all supported Python versions and any other information you want.

2. Ensure that the version of ``cibuildwheel`` used in ``wheels.yml`` is up to date. It may be the case that the latest version of Python is not supported by whichever version is currently in use in ``wheels.yml``.
In this case, ``cibuildwheel`` will not throw an error or warning, it just won't build wheels for that version of Python.

3. Check that the version of ``manylinux`` in the ``cibuildwheel`` section of ``pyproject.toml`` is up to date. This is particularly relevant if any of the errors you encounter are related to specific versions of ``glibc``.

This workflow also builds the source distribution for the package. It is **highly recommended** to download the built sdist from this workflow and test that you can build it against **every** Python version supported by scikit-bio.
Again, this is tedious, but well worth the benefit of being certain things will work before publishing the release.


Tag the release and publish
---------------------------

Assuming that all of the above checks and tests have been performed and are passing, it is time to publish the release.

**Once you hit the ``Publish`` button, all of the automated workflows will run, and they will actually upload the artifacts (documentation, sdist, wheels) to their respective locations (website, PyPI), so ensure that everything is in order before hitting ``Publish``!**

From :repo:`scikit-bio GitHub repo`, click on the releases tab and draft a new release. Use the version number for the tag name (x.y.z) and create the tag against main. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea). This release summary will be the primary information that we point users to when we announce the release.

Hit ``Publish`` when you are ready and keep an eye on the processes to see if any errors arise.


Add wheels to GitHub release page
---------------------------------

Once the automated workflows have run, and our wheels are successfully uploaded to PyPI, you can manually download the wheels from PyPI and add them to the GitHub release page for the release we just published. Click the ``edit release`` button, and then simply add the binaries to the release page.

This isn't an elegant solution, and it may be automated in the future, but for now this is what we're doing.


Post-release cleanup
--------------------

1. Submit and merge a pull request to update the version strings from x.y.z to a.b.c-dev (``skbio.__version__`` should be the only thing needing an update). Update ``CHANGELOG.md`` to include a new section for a.b.c-dev (there won't be any changes to note here yet).

2. Close the release milestone on the GitHub issue tracker if there was one.

3. Send an email to the skbio developers list and anyone else who might be interested (e.g., lab mailing lists). You might include links to the GitHub release page.

4. Tweet about the release from ``@scikit-bio``, including a link to the GitHub release page (for example, https://github.com/scikit-bio/scikit-bio/releases/tag/x.y.z).

5. Beers! :fa:`beer-mug-empty;fa-2x sd-text-success`
