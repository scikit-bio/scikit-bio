Contributing to scikit-bio
==========================

[scikit-bio](http://scikit-bio.org) is an open source software package and we welcome community contributions. You can find the scikit-bio source code on GitHub [here](https://github.com/biocore/scikit-bio).

This document covers what you should do to get started with contributing to scikit-bio. You should read the entire document before contributing code to scikit-bio. This will save time for both you and the scikit-bio developers.

Types of contributions
----------------------

We're interested in many different types of contributions, including feature additions, bug fixes, continuous integration improvements, and documentation/website updates, additions, and fixes.

When considering contributing to scikit-bio, you should begin by posting an issue to the [scikit-bio issue tracker](https://github.com/biocore/scikit-bio/issues). The information that you include in that post will differ based on the type of contribution. Your contribution will also need to be fully tested where applicable (discussed further below).

* For feature additions, please describe why the functionality that you are proposing to add is relevant. For it to be relevant, it should be demonstrably useful to scikit-bio users and it should also fit within the biology/bioinformatics domain. This typically means that a new analytic method is implemented (you should describe why it's useful, ideally including a link to a paper that uses this method), or an existing method is enhanced (e.g., improved performance). We will request benchmark results comparing your method to the pre-existing methods (which would also be required for publication of your method) so pointing to a paper or other document containing benchmark results, or including benchmark results in your issue, will speed up the process. Before contributing a new feature, it's also a good idea to check whether the functionality exists in other Python packages, or if the feature would fit better in another Python package. For example, low-level statistical methods/tests may fit better in a project that is focused on statistics (e.g., [SciPy](http://scipy.org/) or [statsmodels](http://statsmodels.sourceforge.net/)).

* For bug fixes, please provide a detailed description of the bug so other developers can reproduce it. We take bugs in scikit-bio very seriously. Bugs can be related to errors in code, documentation, or tests. Errors in documentation or tests are usually updated in the next scheduled release of scikit-bio. Errors in code that could result in incorrect results or inability to access certain functionality may result in a bug fix release of scikit-bio that is released ahead of schedule.

 You should include the following information in your bug report:

 1. The exact command(s) necessary to reproduce the bug.
 2. A link to all necessary input files for reproducing the bug. These files should only be as large as necessary to create the bug. For example, if you have an input file with 10,000 FASTA-formatted sequences but the error only arises due to one of the sequences, create a new FASTA file with only that sequence, run the command that was giving you problems, and verify that you still get an error. Then post that command and link to the trimmed FASTA file. This is *extremely* useful to other developers and it is likely that if you don't provide this information you'll get a response asking for it. Often this process helps you to better understand the bug as well.

* For documentation additions, you should first post an issue describing what you propose to add, where you'd like to add it in the documentation, and a description of why you think it's an important addition. For documentation improvements and fixes, you should post an issue describing what is currently wrong or missing and how you propose to address it. For more information about building and contributing to scikit-bio's documentation, see our [documentation guide](doc/README.md).

When you post your issue, the scikit-bio developers will respond to let you know if we agree with the addition or change. It's very important that you go through this step to avoid wasting time working on a feature that we are not interested in including in scikit-bio. **This initial discussion with the developers is important because scikit-bio is rapidly changing, including complete re-writes of some of the core objects. If you don't get in touch first you could easily waste time by working on an object or interface that is deprecated.**

Getting started
---------------

### "quick fixes"

Some of our issues are labeled as ``quick fix``. Working on [these issues](https://github.com/biocore/scikit-bio/issues?q=is%3Aopen+is%3Aissue+label%3A%22quick+fix%22) is a good way to get started with contributing to scikit-bio. These are usually small bugs or documentation errors that will only require one or a few lines of code to fix. Getting started by working on one of these issues will allow you to familiarize yourself with our development process before committing to a large amount of work (e.g., adding a new feature to scikit-bio). Please post a comment on the issue if you're interested in working on one of these "quick fixes".

### Joining development

Once you are more comfortable with our development process, you can check out the [``on deck`` label](https://github.com/biocore/scikit-bio/labels/on%20deck) on our issue tracker. These issues represent what our current focus is in the project. As such, they are probably the best place to start if you are looking to join the conversation and contribute code.

Code review
-----------

When you submit code to scikit-bio, it will be reviewed by one or more scikit-bio developers. These reviews are intended to confirm a few points:

* Your code provides relevant changes or additions to scikit-bio ([Types of contributions](#types-of-contributions)).
* Your code adheres to our coding guidelines ([Coding guidelines](#coding-guidelines)).
* Your code is sufficiently well-tested ([Testing guidelines](#testing-guidelines)).
* Your code is sufficiently well-documented ([Documentation guidelines](#documentation-guidelines)).

This process is designed to ensure the quality of scikit-bio and can be a very useful experience for new developers.

Particularly for big changes, if you'd like feedback on your code in the form of a code review as you work, you should request help in the issue that you created and one of the scikit-bio developers will work with you to perform regular code reviews. This can greatly reduce development time (and frustration) so we highly recommend that new developers take advantage of this rather than submitting a pull request with a massive amount of code. That can lead to frustration when the developer thinks they are done but the reviewer requests large amounts of changes, and it also makes it harder to review.

Submitting code to scikit-bio
-----------------------------

scikit-bio is hosted on [GitHub](http://www.github.com), and we use GitHub's [Pull Request](https://help.github.com/articles/using-pull-requests) mechanism for reviewing and accepting submissions. You should work through the following steps to submit code to scikit-bio.

1. Begin by [creating an issue](https://github.com/biocore/scikit-bio/issues) describing your proposed change (see [Types of contributions](#types-of-contributions) for details).

2. [Fork](https://help.github.com/articles/fork-a-repo) the scikit-bio repository on the GitHub website.

3. Clone your forked repository to the system where you'll be developing with ``git clone``. ``cd`` into the ``scikit-bio`` directory that was created by ``git clone``.

4. Ensure that you have the latest version of all files. This is especially important if you cloned a long time ago, but you'll need to do this before submitting changes regardless. You should do this by adding scikit-bio as a remote repository and then pulling from that repository. You'll only need to run the ``git remote`` command the first time you do this:

 ```
 git remote add upstream https://github.com/biocore/scikit-bio.git
 git checkout master
 git pull upstream master
 ```

5. Install scikit-bio for development. See [Setting up a development environment](#setting-up-a-development-environment).

6. Create a new topic branch that you will make your changes in with ``git checkout -b``:

 ```
 git checkout -b my-topic-branch
 ```

 What you name your topic branch is up to you, though we recommend including the issue number in the topic branch, since there is usually already an issue associated with the changes being made in the pull request. For example, if you were addressing issue number 42, you might name your topic branch ``issue-42``.

7. Run ``make test`` to confirm that the tests pass before you make any changes.

8. Make your changes, add them (with ``git add``), and commit them (with ``git commit``). Don't forget to update associated tests and documentation as necessary. Write descriptive commit messages to accompany each commit. We recommend following [NumPy's commit message guidelines](http://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message), including the usage of commit tags (i.e., starting commit messages with acronyms such ``ENH``, ``BUG``, etc.).

9. Please mention your changes in [CHANGELOG.md](CHANGELOG.md). This file informs scikit-bio *users* of changes made in each release, so be sure to describe your changes with this audience in mind. It is especially important to note API additions and changes, particularly if they are backward-incompatible, as well as bug fixes. Be sure to make your updates under the section designated for the latest development version of scikit-bio (this will be at the top of the file). Describe your changes in detail under the most appropriate section heading(s). For example, if your pull request fixes a bug, describe the bug fix under the "Bug fixes" section of [CHANGELOG.md](CHANGELOG.md). Please also include a link to the issue(s) addressed by your changes. See [CHANGELOG.md](CHANGELOG.md) for examples of how we recommend formatting these descriptions.

10. When you're ready to submit your code, ensure that you have the latest version of all files in case some changed while you were working on your edits. You can do this by merging master into your topic branch:

 ```
 git checkout master
 git pull upstream master
 git checkout my-topic-branch
 git merge master
 ```

11. Run ``make test`` to ensure that your changes did not cause anything expected to break.

12. Once the tests pass, you should push your changes to your forked repository on GitHub using:

 ```
 git push origin my-topic-branch
 ```

13. Issue a [pull request](https://help.github.com/articles/using-pull-requests) on the GitHub website to request that we merge your branch's changes into scikit-bio's master branch. Be sure to include a description of your changes in the pull request, as well as any other information that will help the scikit-bio developers involved in reviewing your code. Please include ``fixes #<issue-number>`` in your pull request description or in one of your commit messages so that the corresponding issue will be closed when the pull request is merged (see [here](https://help.github.com/articles/closing-issues-via-commit-messages/) for more details). One of the scikit-bio developers will review your code at this stage. If we request changes (which is very common), *don't issue a new pull request*. You should make changes on your topic branch, and commit and push them to GitHub. Your pull request will update automatically.

Setting up a development environment
------------------------------------

**Note:** scikit-bio must be developed in a Python 3.4 or later environment.

The recommended way to set up a development environment for contributing to scikit-bio is using [Anaconda](https://store.continuum.io/cshop/anaconda/) by Continuum Analytics, with its associated command line utility `conda`. The primary benefit of `conda` over `pip` is that on some operating systems (ie Linux), `pip` installs packages from source. This can take a very long time to install Numpy, scipy, matplotlib, etc. `conda` installs these packages using pre-built binaries, so the installation is much faster. Another benefit of `conda` is that it provides both package and environment management, which removes the necessity of using `virtualenv` separately. Not all packages are available using `conda`, therefore our strategy is to install as many packages as possible using `conda`, then install any remaining packages using `pip`.

1. Install Anaconda

 See [Continuum's site](https://store.continuum.io/cshop/anaconda/) for instructions. [Miniconda](http://conda.pydata.org/docs/install/quick.html) provides a fast way to get conda up and running.

2. Create a new conda environment
 ```
 conda create -n env_name python=3.4 pip
 ```

 Note that `env_name` can be any name desired, for example

 ```
 conda create -n skbio python=3.4 pip
 ```

3. Activate the environment

 This may be slightly different depending on the operating system. Refer to the Continuum site to find instructions for your OS.
 ```
 source activate env_name
 ```

4. Navigate to the scikit-bio directory
 See [the section on submitting code](#submitting-code-to-scikit-bio).
 ```
 cd /path/to/scikit-bio
 ```

5. Install `conda` requirements
 ```
 conda install --file ci/conda_requirements.txt
 ```

6. Install `pip` requirements
 ```
 pip install -r ci/pip_requirements.txt
 ```

7. Install scikit-bio
 ```
 pip install --no-deps -e .
 ```

8. Test the installation
 ```
 make test
 ```

Coding guidelines
-----------------

We adhere to the [PEP 8](http://www.python.org/dev/peps/pep-0008/) Python style guidelines. Please see scikit-bio's [coding guidelines](http://scikit-bio.org/docs/latest/development/coding_guidelines.html) for more details. Before submitting code to scikit-bio, you should read this document carefully and apply the guidelines in your code.

Testing guidelines
------------------

All code that is added to scikit-bio must be unit tested, and the unit test code must be submitted in the same pull request as the library code that you are submitting. We will only merge code that is unit tested and that passes the [continuous integration build](https://github.com/biocore/scikit-bio/blob/master/.travis.yml). This build includes, but is not limited to, the following checks:

- Full unit test suite and doctests execute without errors in supported versions of Python 3.
- C code can be correctly compiled.
- Cython code is correctly generated.
- All tests import functionality from the appropriate minimally deep API.
- Documentation can be built.
- Current code coverage is maintained or improved.
- Code passes ``flake8`` checks.

Running ``make test`` locally during development will include a subset of the full checks performed by Travis-CI.

The scikit-bio coding guidelines describe our [expectations for unit tests](http://scikit-bio.org/docs/latest/development/coding_guidelines.html). You should review the unit test section before working on your test code.

Tests can be executed by running ``make test`` from the base directory of the project or from within a Python or IPython session:

``` python
>>> from skbio.test import pytestrunner
>>> pytestrunner()
# full test suite is executed
```

Documentation guidelines
------------------------

We strive to keep scikit-bio well-documented, particularly its public-facing API. See our [documentation guide](doc/README.md) for more details.

Getting help with git
---------------------

If you're new to ``git``, you'll probably find [gitref.org](http://gitref.org/) helpful.
