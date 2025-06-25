Contribute to scikit-bio
========================

**Scikit-bio** is a community-driven open-source software project, and we warmly welcome your contributions!

We are interested in many types of contributions, including feature additions, bug fixes, continuous integration improvements, and documentation/website updates, additions, and fixes. Whether you are a researcher, educator, or developer; whether your interest lies in biology, mathematics, statistics, or computer science; your input is invaluable. You can help making scikit-bio a better software package for our entire community.

This document covers the information you may need to get started with contributing to scikit-bio. In addition, for a broader perspective, we recommend the inspiring guide: `How to Contribute to Open Source <https://opensource.guide/how-to-contribute/>`_.

Visit our GitHub repository: :repo:`scikit-bio/scikit-bio <>` for the source code of scikit-bio. You will need a GitHub account to interact with the scikit-bio codebase and community.

Valuable contributions can be made without or with minimum amount of coding. We detail various ways you can contribute below:

- `Ask a question`_ | `Report an error`_ | `Suggest a new feature`_ | `Fix a typo`_

Contributing code to scikit-bio is a rigorous and rewarding process. We have prepared the following step-by-step guidelines:

- `Before coding`_ | `Set up a workspace`_ | `Write code`_ | `Test code`_ | `Document code`_ | `Style code`_ | `Submit code`_ | `Review code`_

.. .. contents::
..    :depth: 1
..    :local:
..    :backlinks: none

In addition, there are separate documents covering advanced topics:

.. toctree::
   :maxdepth: 1

   devdoc/code_guide
   devdoc/doc_guide
   devdoc/new_module
   devdoc/review
   devdoc/release


Ask a question
--------------

Your inquiry matters! By asking questions in the scikit-bio :repo:`issue tracker <issues>` or :repo:`discussion board <discussions>`, you are not only giving us (and the community) the chance to help, but also let us assess the needs of users like you. Before asking a question, take a moment to search existing threads to see if there are any relevant ones. We also keep an eye on broader community forums such as Stack Overflow and BioStars for questions related to our scope.


Report an error
---------------

The scikit-bio team is proud of our high-quality, well-tested codebase. That being said, no software is immune to errors, which may arise from bugs, overlooked edge cases, or confusions in documentation. In any situation, we would appreciate it if you can report the error you encountered to us.

You may :repo:`open an issue <issues/new/choose>` to report the error. Please provide a detailed description of the error such that the developers can reproduce it. Specifically, you may include the following information in the report:

1. The exact **command(s)** necessary to reproduce the error.

2. The input **file(s)** necessary for reproducing the error. You may either attach the file in the issue (by dragging & dropping) if it is small, or provide a link to it otherwise. The file should only be as large as necessary to reproduce the error.

.. note:: For example, if you have a FASTA file with 10,000 sequences but the error only arises due to one of the sequences, create a new FASTA file with only that sequence, run the command that was giving you problems, and verify that you still get an error. Then post that command and link to the trimmed FASTA file.

This is *extremely* useful to the developers, and it is likely that if you don't provide this information you'll get a response asking for it. Often this process helps you to better understand the error as well.

We take error reports very seriously. Once confirmed that they should be fixed, we will update the code to fix them as soon as we can, and ship the update in the next scheduled release of scikit-bio. If the error could result in incorrect results or inability to access certain functionality, we may release a bug-fix version of scikit-bio ahead of the schedule.


Suggest a new feature
---------------------

We are always looking for new ideas to enhance scikit-bio's capabilities, especially from users with unique research interests. If you believe there is an analysis or feature that could extend scikit-bio's current offerings, we warmly invite you to share your suggestions with us.

Please describe why the functionality that you are suggesting is relevant. For it to be relevant, it should be demonstrably useful to scikit-bio users and it should also fit within the biology/bioinformatics domain. This typically means that a new analytic method is implemented (you should describe why it's useful, ideally including a link to a paper that uses this method), or an existing method is enhanced (e.g., improved performance).

If the scope of the suggested method overlaps with any pre-existing methods in scikit-bio, we may request benchmark results comparing your method to the pre-existing ones (which would also be required for publication of your method) so pointing to a paper or other document containing benchmark results, or including benchmark results in your issue, will help.

Before suggesting a new feature, it is also a good idea to check whether the functionality exists in other Python packages, or if the feature would fit better in another Python package. For example, low-level statistical methods/tests may fit better in a project that is focused on statistics (e.g., `SciPy <https://scipy.org/>`_ or `statsmodels <https://www.statsmodels.org/>`_).

If your proposal represents a significant research direction or requires a substantial suite of methods, we encourage you to consider establishing a formal academic or industrial collaboration with the scikit-bio team. For more details on this process, please refer to the :ref:`about:Collaboration` section.


Fix a typo
----------

If you spot small errors such as typos, redundant spaces, broken links, missing citations etc. in the scikit-bio code or documentation, and want to give it a quick fix, you may follow the procedures detailed below. All procedures will take place in the web browser, and don't involve creating anything in your local computer.

.. warning:: This approach should not be applied to anything larger than small errors. For the latter, please read `Before coding`_.

1. Locate the specific code file that needs to be fixed in the GitHub repository. If you are reading the documentation, you can click the `[source] <about:blank>`__ link next to the header to locate the corresponding code file.

2. In the top-right corner of the code viewer there is an Edit (:octicon:`pencil`) button, with a prompt "Fork this repository and edit the file". Click it. Then click the button :bdg-success:`Fork this repository`. This will open GitHub's `online file editor <https://docs.github.com/en/repositories/working-with-files/managing-files/editing-files>`_.

3. **Edit the code**.

4. When done, click :bdg-success:`Commit changes...`. Then enter a commit message to describe what you did, like "Fixed a typo in the documentation of skbio.module.function". Then click :bdg-success:`Propose changes`.

5. You will be able to review the changes you made and compare with the original code :octicon:`git-compare`. If everything looks good to you, click :bdg-success:`Create pull request`. Then enter a title and description that you think are informative to the scikit-bio maintainers. The **title** may or may not be the same as the commit message. In the **description**, you will need to answer a few questions by typing an ``x`` in the relevant checkboxes. You may also explain why the original code should be replaced by yours. Finally, click :bdg-success:`Create pull request`.

6. This will create a `pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`__ :octicon:`git-pull-request`, i.e., the changes you made to the scikit-bio repository. A scikit-bio maintainer will review your pull request, and run necessary tests to make sure it is sound. You may be asked to clarify or to make modifications to your code. Please work with the maintainer by replying in the pull request.

7. If the maintainer believes that your code is good to go, they will `merge <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request>`_ it into the scikit-bio codebase. Once merged, the pull request webpage will have a purple notice :octicon:`git-merge`, saying: "Pull request successfully merged and closed".

8. At this point, your contribution is completed 🎉. You may optionally click :bdg-light:`Delete branch` to clean up your workspace. Then, you can move on, and enjoy the improved scikit-bio!


Before coding
-------------

We sincerely value your willingness to contribute code to scikit-bio (beyond reporting issues or correcting typos). This process can be intensive, particularly for those new to software engineering. The following sections detail the steps for contributing code to scikit-bio. Please review them carefully.

Discuss your plan
^^^^^^^^^^^^^^^^^

When considering contributing code to scikit-bio, you should begin by posting an issue to the :repo:`scikit-bio issue tracker <issues>`. The information that you include in that post will differ based on the type of contribution. The two types of issues discussed in `Report an error`_ and `Suggest a new feature`_ can be a good start of the discussion.

The scikit-bio developers will respond to let you know if we agree with the addition or change. It's very important that you go through this step to avoid spending time working on a feature that we are not interested in including in scikit-bio.

Take existing tasks
^^^^^^^^^^^^^^^^^^^

Alternatively, if you're looking to contribute where help is needed, you can explore the following types of issues:

- **Quick fix**: Some of our issues are labeled as ``quick fix``. Working on :repo:`these issues <issues?q=is%3Aopen+is%3Aissue+label%3A%22quick+fix%22>` is a good way to get started with contributing to scikit-bio. These are usually small bugs or documentation errors that will only require one or a few lines of code to fix. Getting started by working on one of these issues will allow you to familiarize yourself with our development process before committing to a large amount of work (e.g., adding a new feature to scikit-bio). Please post a comment on the issue if you're interested in working on one of these "quick fixes".

- **On deck**: Once you are more comfortable with our development process, you can check out the ``on deck`` :repo:`label <labels/on%20deck>` on our issue tracker. These issues represent what our current focus is in the project. As such, they are probably the best place to start if you are looking to join the conversation and contribute code.


Set up a workspace
------------------

To start contributing code to scikit-bio, you'll need to prepare a local development environment. This section guides you through the process step-by-step.

1. `Fork <https://help.github.com/articles/fork-a-repo>`_ the scikit-bio repository on the GitHub website. This will create a copy of the repository under your account, and you can access it using the URL: ``https://github.com/urname/scikit-bio/`` (``urname`` is your GitHub account).

2. `Clone <https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository>`_ the forked repository to your local computer. Specifically, in your forked repository, you may click the :bdg-success:`<> Code` button, make sure that you are under the "Local" - "SSH" tab, and copy the URL to the clipboard. It should typically be: ``git@github.com:urname/scikit-bio.git``.

.. note:: If the "SSH" tab is not available, it could mean that you have not set up an SSH key for your GitHub account. You may follow the instructions on `Connecting to GitHub with SSH <https://docs.github.com/en/authentication/connecting-to-github-with-ssh>`_ to set it up.

Then, open "Terminal" :octicon:`terminal` (or anything similar) in your local computer and navigate to a directory where you want to place the workspace, and execute::

    git clone git@github.com:urname/scikit-bio.git

.. note::

   If this is the first time you use ``git``, you may follow the `Set up Git <https://docs.github.com/en/get-started/getting-started-with-git/set-up-git>`_ guidelines to install Git and set your user name and email address.

   This tutorial assumes that you will use the classical ``git`` to create a local development environment. If you prefer other methods such as `GitHub CLI <https://cli.github.com/>`_ or `Codespace <https://github.com/features/codespaces>`_, please follow corresponding instructions.

This will create a directory ``scikit-bio`` containing all files in the repository. Enter the directory::

    cd scikit-bio

Add the official scikit-bio repo as the **upstream** of your fork::

    git remote add upstream https://github.com/scikit-bio/scikit-bio.git

3. Create a development environment with necessary dependencies. This is typically done using `Conda <https://conda.io/>`_ (or `Mamba <https://mamba.readthedocs.io/en/latest/index.html>`_, in which case the command ``conda`` in the following code should be replaced with ``mamba``).

.. note::

   If you do not have Conda (or Mamba) in your computer, you may install one of the distributions such as `Miniconda <https://conda.io/miniconda.html>`_, `Miniforge <https://github.com/conda-forge/miniforge>`_ or `Anaconda <https://www.anaconda.com/download/>`_.

   We recommend Conda over other approaches such as ``pip``, ``pyenv``, and ``virtualenv``. However, you are not blocked from using them in necessary situations.

Execute the following command (``skbio-dev`` can be any name you like)::

    conda create -n skbio-dev -c conda-forge --file ci/conda_requirements.txt --file ci/requirements.test.txt --file ci/requirements.lint.txt --file ci/requirements.doc.txt

When done, activate the environment::

    conda activate skbio-dev

.. note:: This may be slightly different depending on the operating system. Refer to the `Conda documentation <https://docs.conda.io/>`_ to find instructions for your OS.

4. Install scikit-bio from source code::

    pip install --no-deps -e .

This will install scikit-bio to the current conda environment. After this, you can use scikit-bio like a normal user (e.g., you can do ``import skbio`` in Python code). When you edit the code in the this directory, the changes will be immediately reflected as you use the software.

5. Test the installation::

    make test

This will run all unit tests implemented in the scikit-bio codebase to check if the corresponding functionality works correctly. The output should only indicate passes and warnings, but no failures.

6. Activate pre-commit hooks::

    pre-commit install

This will enable a set of tools that will automatically execute every time you commit changes to ensure code quality.


Write code
----------

Before you start writing code, you may discuss with the scikit-bio team to make sure that your intended contribution is relevant (see `Before coding`_ above). Next, you may work through the following steps to start coding.

1. Update your main branch such that it has the latest version of all files. This is especially important if you cloned a long time ago::

    git checkout main
    git pull upstream main

Optionally, you may do the following to keep your forked repository's main branch up-to-date as well::

    git push origin main

2. Create a new `branch <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches>`__ that you will make changes in::

    git checkout -b mywork

``mywork`` is the name of the branch. What you name the branch is up to you, though we recommend including the issue number, if there is a relevant one (see above). For example, if you were addressing issue #42, you might name your branch ``issue-42``.

.. warning:: It is not recommended that you directly code in your ``main`` branch.

3. Run ``make test`` to confirm that the tests pass before you make any changes.

4. **Edit the code** using any code editor of your favor. Now it is the time to bring your creativity to scikit-bio!

Scikit-bio's :doc:`coding guidelines <devdoc/code_guide>` provide more details on how to write high-quality code. It is recommended that you read this document carefully and apply the guidelines in your code.


Test code
---------

Testing your code is an essential step to ensure its quality and functionality. For scikit-bio, we emphasize the importance of comprehensive testing to guarantee that your contribution works as expected. You might find testing tedious (it often costs more time than implementing the algorithm!), but it is a valuable step that will significantly enhance the robustness of your code, and will help fostering a culture of reliability and trust within our community. This section is to guide you through the testing process, providing you with the tools and knowledge to perform effective tests for scikit-bio.

Functional test
^^^^^^^^^^^^^^^

You will want to test your code like a *user* would do: import the function, execute it on some input data, and examine whether the results are correct. You may install additional software such as `Jupyter <https://jupyter.org/>`_ in the same conda environment to make the testing process convenient and pleasant. There is no need to reinstall the modified scikit-bio package in a separate environment in order to test. As soon as you edit the code, the changes will immediately reflect when you use the code.

For example, you wrote a function ``gc_content``, which calculates the fraction of ``G`` and ``C`` in a nucleotide sequence::

    def gc_content(seq):
        return sum(1 for x in seq if x in "GC") / len(seq)

You added this function to the code file ``skbio/sequence/nucl.py``. It will be available for use in any Python code launched in the same conda environment::

    >>> from skbio.sequence.nucl import gc_content

Now test the function on some data. For example, you would expect that ``gc_content("ACGT")`` and ``gc_content("GGATCCGC")`` return 0.5 and 0.75, respectively. Is that the case?

We highly recommend that you use real-world biological data in addition to small, dummy data for testing. This will let you evaluate the robustness and scalability of your code in real-world applications. For example, DNA sequences retrieved from public databases may contain lowercase letters, and you will find that the ``gc_content`` function cannot handle them properly. For another instance, a FASTQ file with ten million sequences (which is common) may cost a function forever to process, in which case you should consider optimization.

Unit test
^^^^^^^^^

`Unit testing <https://en.wikipedia.org/wiki/Unit_testing>`_ involves testing the smallest units of your code, such as classes, functions, and methods, to ensure they function correctly in isolation. It is a fundamental best practice in software engineering, but is often overlooked by beginners. Unit testing is made easier by writing test code alongside the algorithm code. Both types of code are integrated into the scikit-bio codebase. This test code is then regularly executed whenever changes are made to ensure that the intended behavior remains consistent over time.

For example, the test code for the ``gc_content`` may live in ``skbio/sequence/tests/test_nucl.py``, under class ``TestNucl``, as a method ``test_gc_content``. It may read like::

    def test_gc_content(self):
        self.assertEqual(gc_content("ACGT"), 0.5)
        self.assertEqual(gc_content("GGATCCGC"), 0.75)
        ...

You can run this test with::

    python skbio/sequence/tests/test_nucl.py TestNucl.test_gc_content

The screen output will tell you whether the test passed, and if not, what went wrong. This information will help you debug your code.

Ideally, every line of the code should be covered by unit test(s). For example, if your function has an ``if`` statement, both ``True`` and ``False`` situations should be tested.

It is a good practice to test all types of cases you can think of, including normal cases and `edge cases <https://en.wikipedia.org/wiki/Edge_case>`_. For example, an empty sequence (``""``) will cause the ``gc_content`` function to crash, because zero cannot serve as a denominator in an equation. Having edge cases like this will help you to identify limitations of your code and think whether you should implement special handling to avoid problems.

You should also test whether the changed code fits into scikit-bio without causing problems in the other parts of the codebase. There is a convenient command to run all unit tests implemented in scikit-bio::

    make test

Alternatively, you may run all unit tests in a Python session (including Jupyter)::

    >>> from skbio.test import pytestrunner
    >>> pytestrunner()

Code coverage
^^^^^^^^^^^^^

`Code coverage <https://en.wikipedia.org/wiki/Code_coverage>`_ refers to the percentage of source code lines covered by unit tests. It is an assessment of the quality of a software project. In scikit-bio, code coverage can be calculated using the following command::

    coverage run -m skbio.test && coverage report

This will report the coverage of each code file and the entire codebase. If the coverage decreased for the file you edited, you may have missed some anticipated unit tests. You can create a detailed HTML report with::

    coverage html

Then open ``htmlcov/index.html`` in a web browser, navigate to the page for the relevant code file, and check which lines of your code are not covered by unit tests. Work on them to bring back coverage.

Please read the :ref:`devdoc/code_guide:How should I test my code?` section of the coding guidelines to learn more about unit testing.


Document code
-------------

`Documentation <https://en.wikipedia.org/wiki/Software_documentation>`_ is a vital part of software engineering, especially for projects like scikit-bio, which involve many contributors and are designed to endure over time. Documentation helps everyone -- users and developers -- get on the same page. It also helps you, as even seasoned developers can lose track of their own coding logic over time. Also remember that scikit-bio brings together people from various fields, and nobody is expected to have the same level of understanding across all disciplines. Therefore, documenting your code with the broader audience in mind is important. This section will cover the basics of documenting your code in a manner that benefits the scikit-bio community at large.

Scikit-bio's :doc:`documentation guidelines <devdoc/doc_guide>` provide more details on how to write effective documentation.

Comments
^^^^^^^^

`Comments <https://en.wikipedia.org/wiki/Comment_(computer_programming)>`_ in the source code explain the rationale to fellow developers. Please make comments frequently in your code, especially where the code itself is not that intuitive. For example::

    # Perform eigendecomposition on the covariance matrix of data.
    # w and v are eigenvalues and eigenvectors, respectively.
    # eigh is used in favor of eig to avoid returning complex numbers due to
    # matrix asymmetry caused by floating point errors.
    # Discussed in: https://stackoverflow.com/questions/12345678/
    w, v = np.linalg.eig(np.cov(X))

Please read the :ref:`devdoc/code_guide:How should I write comments?` section of the coding guidelines to learn more about writing comments.

Docstrings
^^^^^^^^^^

`Docstrings <https://en.wikipedia.org/wiki/Docstring>`_ are structured text blocks associated with each unit of the code that detail the usage of the code. Docstrings will be rendered to the software documentation. That is, *users* (not just developers) will read them. Therefore, docstrings are critical if you want your code to be used, and in the correct way.

Below is a very simple example for the ``gc_content``. The lines between the triple double quotes (``"""``) are the docstring::

    def gc_content(seq):
        """Calculate the GC content of a nucleotide sequence.

        Parameters
        ----------
        seq : str
            Input sequence.

        Returns
        -------
        float
            Fraction of G and C.

        """
        return sum(1 for x in seq if x in "GC") / len(seq)

As shown, the docstring explains the purpose, the parameter(s), and the return value(s) of the function. In more complicated cases, the docstring should also include example usage, potential errors, related functions, mathematics behind the algorithm, references to webpages or literature, etc. Every public-facing component should have a docstring.

Please read the :ref:`devdoc/doc_guide:Docstring style` section of the documentation guidelines to learn more about writing docstrings.

Doctests
^^^^^^^^

You may consider adding **example usages** of your code to its docstring. For example::

    def gc_content(seq):
        """Calculate the GC content of a nucleotide sequence.
        ...

        Examples
        --------
        >>> from skbio.sequence.nucl import gc_content
        >>> gc_content("ACGT")
        0.5

        """
        ...

The example code and its output must match. This is ensured by `doctest <https://docs.python.org/3/library/doctest.html>`_. When you run ``make test`` (see above), doctests are automatically executed as part of the test suite. You may fix any issues according to the screen output.

HTML rendering
^^^^^^^^^^^^^^

After completing docstrings, you will want to check how they look like when rendered to the documentation webpages. You may build the entire HTML documentation package locally with::

    make doc

The built documentation will be at ``doc/build/html/index.html``, and can be examined using your web browser. If errors arise during the building process, or the rendered webpages don't look as anticipated, you should address the issues accordingly.

Changelog
^^^^^^^^^

Please mention your changes in :repo:`CHANGELOG.md <blob/main/CHANGELOG.md>`. This file informs scikit-bio *users* of changes made in each release, so be sure to describe your changes with this audience in mind. It is especially important to note API additions and changes, particularly if they are backward-incompatible, as well as bug fixes. Be sure to make your updates under the section designated for the latest development version of scikit-bio (this will be at the top of the file). Describe your changes in detail under the most appropriate section heading(s). For example, if your pull request fixes a bug, describe the bug fix under the "Bug fixes" section of the changelog. Please also include a link to the issue(s) addressed by your changes.


Style code
----------

`Code style <https://en.wikipedia.org/wiki/Programming_style>`_ is a set of rules for formatting and structuring code in a particular software project. Although violating these rules won't cause errors in executing the code, adhering to them ensures that the codebase remains uniform and professional, and facilitates team collaboration.

Scikit-bio utilizes the `Ruff <https://docs.astral.sh/ruff/>`_ program for autoformatting and linting to ensure code consistency and quality. The rules are specified in :repo:`pyproject.toml <blob/main/pyproject.toml>`. Basically, we largely adopted the `Black <https://black.readthedocs.io/>`_ code style.

When you `set up the development environment <#set-up-a-workspace>`_, Ruff was already installed and integrated into a `pre-commit hook <https://github.com/astral-sh/ruff-pre-commit>`__. This means that Ruff will automatically format and lint your code every time you commit changes (see `Submit code`_ below). Therefore *you are not required to take any explicit action*. However, you can still manually run Ruff to check and fix issues in specific code files you have worked on::

    ruff check --fix mycode.py

If Ruff identifies any errors that cannot be automatically fixed, you will need to manually fix them based on Ruff's feedback. When done, let Ruff reformat your code::

    ruff format mycode.py

You will notice the improvement in your code's appearance before and after using Ruff. While it is always beneficial to strive for professional-looking code from the start, the necessity for perfection has lessened with the advent of tools like Ruff.


Submit code
-----------

Having completed, tested, and documented your code, you may now believe it deserves a place in scikit-bio to benefit the community. This section outlines the steps to submit your code to the official scikit-bio repository.

1. Add any new code file(s) you created to the git repository::

    git add path/to/mycode.py

Alternatively, if you have multiple new files, you can add them all at once::

    git add .

2. `Commit <https://github.com/git-guides/git-commit>`__ the changes (this is like "saving" your code in the current branch)::

    git commit -am "describe what I did"

Here, "describe what I did" is the placeholder of a *commit message*. You should write a meaningful commit message to describe what you did. We recommend following `NumPy's commit message guidelines <https://numpy.org/doc/stable/dev/development_workflow.html#writing-the-commit-message>`_, including the usage of commit tags (i.e., starting commit messages with acronyms such ``ENH``, ``BUG``, etc.).

The `commit` command will trigger the `pre-commit hook <https://pre-commit.com/>`__, which automatically runs Ruff to check and fix any code style problems (see `Style code`_ above). If there are any errors flagged by Ruff, you will need to resolve them and commit again.

3. Merge the latest code from the official scikit-bio repository to your local branch::

    git fetch upstream
    git merge upstream/main

This step is important as it ensures that your code doesn't conflict with any recent updates in the official repository. This could happen as there are other developers simultaneously working on the project.

If there are conflicts, you will need to `resolve the conflicts <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/addressing-merge-conflicts/resolving-a-merge-conflict-on-github>`__ by editing the affected files. When done, run ``git add`` on those files, then run ``git commit`` with a relevant commit message (such as "resolved merge conflicts").

4. Run ``make test`` for the last time to ensure that your changes don't cause anything to break.

5. Once the tests pass, you should push your changes to your forked repository on GitHub using::

    git push origin mywork

6. Navigate to the GitHub website, and create a `pull request <https://help.github.com/articles/using-pull-requests>`__ from your ``mywork`` branch to the ``main`` branch of the official scikit-bio repository. Usually, GitHub will prompt you to do so, and you may click the :bdg-success:`Compare & pull request` button to initiate this process. If not, you can invoke a :bdg-success:`New pull request` under the ":octicon:`git-pull-request` pull request" tab.

7. Enter a meaningful title and a description of your code in the pull request. You may mention the issue you attempt to address in the description, such as "Resolves #42". This will `link your pull request to the issue <https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue>`__. You will also need to answer a few questions by either typing an ``x`` in the checkboxes that apply to your code or leaving them empty otherwise. These questions can be found in :repo:`PULL_REQUEST_TEMPLATE.md <blob/main/.github/PULL_REQUEST_TEMPLATE.md>`. When done, click :bdg-success:`Create pull request`.


Review code
-----------

Your `pull request will be reviewed <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/reviewing-changes-in-pull-requests/about-pull-request-reviews>`__ by one or more maintainers of scikit-bio. These reviews are intended to confirm a few points:

- Your code provides relevant changes or additions to scikit-bio.
- Your code adheres to our coding guidelines.
- Your code is sufficiently well-tested.
- Your code is sufficiently well-documented.

This process is designed to ensure the quality of scikit-bio and can be a very useful experience for new developers.

Typically, the reviewer will launch some automatic checks on your code. These checks are defined in :repo:`ci.yml <blob/main/.github/workflows/ci.yml>`. They involve:

- Full unit test suite and doctests execute without errors in all supported software and hardware environments.
- C code can be correctly compiled.
- Cython code is correctly generated.
- Documentation can be built.
- Code coverage is maintained or improved.
- Code passes linting.

The checks may take several to a few dozen minutes. If some check(s) fail, you may click "Details" in these checks to view the error messages, and fix the issues accordingly.

Meanwhile, the reviewer will comment on your code inline and/or below. They may request changes (which is very common). Please work with the reviewer to improve your code.

You should revise your code in your local branch. When completed, commit and push your code again (steps 1-5 of `Submit code`_). This will automatically update your pull request and restart the checks. *Don't issue a new pull request*.

.. note:: Particularly for big changes, if you'd like feedback on your code in the form of a code review as you work, you should request help in the issue that you created and one of the scikit-bio maintainers will work with you to perform regular code reviews. This can greatly reduce development time. We highly recommend that new developers take advantage of this rather than submitting a pull request with a massive amount of code. That can lead to frustration when the developer thinks they are done but the reviewer requests large amounts of changes, and it also makes it harder to review.

Please read :doc:`devdoc/review` for more details on how pull requests should be reviewed in the scikit-bio project.

After your code has been improved and the reviewer has approved it, they will `merge your pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/incorporating-changes-from-a-pull-request/merging-a-pull-request>`__ into the ``main`` branch of the official scikit-bio repository. This will be indicated by a note: :octicon:`git-merge` "Pull request successfully merged and closed".

Congratulations! Your code is now an integral part of scikit-bio, and will benefit the broader community. You have successfully completed your contribution, and we extend our appreciation to you! 🎉🎉🎉"


Benchmarks
----------

Scikit-bio utilizes `ASV <https://github.com/airspeed-velocity/asv>`_ to run benchmarks on a selection of its functions. Benchmarks help the development team prevent performance regression, the unintentional loss of performance which may come from modified code. Scikit-bio's benchmarks are run against every release of scikit-bio, starting with version ``0.6.1``. We welcome the addition of new benchmarks or the expansion of existing benchmarks to further protect against performance regression. If you are interested in contributing to our benchmarking system please see our `benchmark repository <https://github.com/scikit-bio/scikit-bio-benchmarks>`_.