# Releasing a new version of scikit-bio

## Introduction

This guide explains how to release a new version of scikit-bio. To illustrate examples of commands you might run, let's assume that the current version is x.y.y-dev and we want to release version x.y.z.

**Note:** The following commands assume you are in the top-level directory of the scikit-bio repository unless otherwise noted. They also assume that you have [Miniconda3](http://conda.pydata.org/miniconda.html) installed. It is important that you use Miniconda3, and not Miniconda2, because the root `conda` environment needs Python 3 for the build steps below.

## Prep the release

1. Ensure the Travis build is passing against master.

2. Update the version strings (x.y.y-dev) to the new version (x.y.z). This will include `__version__` defined in ``skbio/__init__.py``, as well as any `@experimental/@stable/@deprecated` [API stability decorators](http://scikit-bio.org/docs/latest/user/api_stability.html) with `as_of='x.y.y-dev'`. ``grep`` for the current version string to find all occurrences:

        grep -r 'x\.y\.y-dev' .

3. Remove any deprecated functionality that was scheduled for removal on or before this release. When removing deprecated functionality, make sure the functionality has been in a deprecated state for the appropriate number of releases described in the [API stability docs](http://scikit-bio.org/docs/latest/user/api_stability.html). If there is functionality that shouldn't be removed yet, bump the `until` version to a future version. To find all deprecated functionality, search for `@deprecated` decorators:

        grep -r '@deprecated' .

    Note any deprecated functionality that was removed in the `### Miscellaneous` section of `CHANGELOG.md`.

4. Update ``CHANGELOG.md`` to include descriptions of all **user-facing** changes that made it into this release. Be sure to update the heading to include the new version (x.y.z) and the date of the release. Use the existing structure in the file as a template/guide.

5. Submit a pull request and merge when Travis-CI tests are passing.

## Build website docs

You will need to **fully install** the latest master branch of scikit-bio (including built extensions) and build the docs from this version. **Make sure the version of scikit-bio that is imported by ``import skbio`` is the correct one before building the docs.**

1. Build the documentation locally:

        make -C doc clean html

2. Switch to the ``gh-pages`` branch of the repository.

3. Remove ``docs/latest``:

        git rm -rf docs/latest

4. Copy over the built documentation to ``docs/x.y.z`` and ``docs/latest``:

        cp -r doc/build/html docs/latest
        cp -r doc/build/html docs/x.y.z

5. Add a new list item to ``index.html`` to link to ``docs/x.y.z/index.html``.

6. Port content from ``README.md`` to ``index.html`` if there are any changes that need to be included on the front page.

7. Test out your changes by opening the site locally in a browser. Be sure to check the error console for any errors.

8. Submit a pull request with the website updates and merge. **Note:** Once merged, the live website is updated, so be sure to poke through the live site to make sure things are rendered correctly with the right version strings.

## Tag the release

From the [scikit-bio GitHub page](https://github.com/biocore/scikit-bio), click on the releases tab and draft a new release. Use the version number for the tag name (x.y.z) and create the tag against master. Fill in a release title that is consistent with previous release titles and add a summary of the release (linking to ``CHANGELOG.md`` is a good idea). This release summary will be the primary information that we point users to when we announce the release.

Once the release is created on GitHub, it's a good idea to test out the release tarball before publishing to PyPI:

1. Create a new `conda` environment for testing (fill in a name for `<environment>`):

        conda create -n <environment> python=3.5 numpy
        source activate <environment>

2. Install the release tarball from GitHub and run the tests:

        pip install https://github.com/biocore/scikit-bio/archive/x.y.z.tar.gz
        python -m skbio.test

## Publish the release

Assuming the GitHub release tarball correctly installs and passes its tests, you're ready to create the source distribution (``sdist``) that will be published to PyPI. It is important to test the source distribution because it is created in an entirely different way than the release tarball on GitHub. Thus, there is the danger of having two different release tarballs: the one created on GitHub and the one uploaded to PyPI.

1. Download the release tarball from GitHub, extract it, and ``cd`` into the top-level directory.

2. Build a source distribution:

        python setup.py sdist

3. Create and activate a new `conda` environment, and test the `sdist`:

        pip install dist/scikit-bio-x.y.z.tar.gz
        cd  # cd somewhere outside the extracted scikit-bio directory
        python -m skbio.test

4. If everything goes well, it is finally time to push the release to PyPI:

        python setup.py sdist upload

    You must have the proper login credentials to add a release to PyPI. Currently [@gregcaporaso](https://github.com/gregcaporaso) has these, but they can be shared with other release managers.

5. Once the release is available on PyPI, do a final round of testing. Create a new `conda` environment and run:

        pip install scikit-bio
        cd  # cd somewhere outside the extracted scikit-bio directory
        python -m skbio.test

    If this succeeds, the PyPI release appears to be a success. Make sure the installed version is the correct one.

6. Next, we'll prepare and post the release to [anaconda.org](http://www.anaconda.org).

    You'll need to have ``conda-build`` and ``anaconda-client`` installed to perform these steps. Both can be conda-installed. First, log into anaconda with your anaconda username using the following command. You should have access to push to the ``biocore`` anaconda account through your account (if you don't, get in touch with [@gregcaporaso](https://github.com/gregcaporaso) who is the owner of that account).

        anaconda login

    Due to its C extensions, releasing scikit-bio packages for different platforms will require you to perform the following steps on each of those platforms. For example, an ``osx-64`` package will need to be built on OS X, and a ``linux-64`` package will need to be built on 64-bit Linux. These steps will be the same on all platforms, so you should repeat them for every platform you want to release for.

        conda skeleton pypi scikit-bio
        conda build scikit-bio --python 3.4
        conda build scikit-bio --python 3.5

    **Note:** When building 64-bit Linux packages, it is recommended that you use conda-forge's `linux-anvil` Docker image. This ensures a consistent Linux build environment that has an old enough version of `libc` to be compatible on most Linux systems. To start up a `linux-anvil` Docker container:

        docker run -i -t condaforge/linux-anvil
        # Now you should be in the linux-anvil environment
        sed -i '/conda-forge/d' ~/.condarc
        # Run the build commands from above

    At this stage you have built Python 3.4 and 3.5 packages. The absolute path to the packages will be provided as output from each ``conda build`` commands. You should now create conda environments for each, and run the tests as described above. You can install these local packages as follows:

        conda install --use-local scikit-bio

    If the tests pass, you're ready to upload.

        anaconda upload -u biocore <package-filepath>

    ``<package-filepath>`` should be replaced with the path to the package that was was created above. Repeat this for each package you created (here, the Python 3.4 and 3.5 packages).

    After uploading, you should create new environments for every package you uploaded, install scikit-bio from each package, and re-run the tests. You can install the packages you uploaded as follows:

        conda install -c https://conda.anaconda.org/biocore scikit-bio

## Post-release cleanup

1. Submit and merge a pull request to update the version strings from x.y.z to x.y.z-dev (`skbio.__version__` should be the only thing needing an update). Update ``CHANGELOG.md`` to include a new section for x.y.z-dev (there won't be any changes to note here yet).

2. Close the release milestone on the GitHub issue tracker if there was one.

3. Send an email to the skbio developers list and anyone else who might be interested (e.g., lab mailing lists). You might include links to the GitHub release page.

4. Tweet about the release from `@scikit-bio`, including a link to the GitHub release page (for example, https://github.com/biocore/scikit-bio/releases/tag/x.y.z). Post a similar message to [scikit-bio's Gitter](https://gitter.im/biocore/scikit-bio).

5. :beers:
