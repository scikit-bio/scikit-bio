|license| |build| |coverage| |bench| |release| |pypi| |conda| |gitter|

.. image:: logos/logo.svg
   :width: 600 px
   :target: https://scikit.bio
   :alt: scikit-bio logo

*scikit-bio is an open-source, BSD-licensed Python 3 package providing data structures, algorithms and educational resources for bioinformatics.*


🌟 Announcing the revitalization of scikit-bio! 🌟
--------------------------------------------------

We are excited to announce the resurgence of **scikit-bio** (`<https://scikit.bio>`_)! With a re-assembled developer team and funding from the DOE, we're back with renewed vigor in 2024!

Our vision is to expand scikit-bio into a more robust and versatile library, catering to the ever-growing demands of multi-omic data analysis. This resurgence marks a new chapter in which we will focus on:

- Streamlining the analysis of diverse, massive omic data, emphasizing efficiency and versatility.
- Integrating advanced techniques for multi-omic integration to unravel the complex interplay between biological systems and environments.
- Implementing methods for modeling and annotating biological features utilizing community ecology and phylogenetics.

We invite the scientific community to join us in shaping the future of scikit-bio. Your expertise, feedback, and contributions will be the driving force behind this exciting phase.

Stay `tuned for updates <https://github.com/scikit-bio/scikit-bio/discussions/categories/announcements>`_, and let's innovate together for a deeper understanding of bio-complexities!


Get involved with scikit-bio
----------------------------

Your questions, ideas, and contributions matter!

Join our community on `GitHub Discussions <https://github.com/scikit-bio/scikit-bio/discussions>`_: This is your go-to place for asking questions, sharing insights, and participating in discussions about scikit-bio. Engage with both the developers and fellow users here.

Report issues and bugs: If you encounter specific problems when using scikit-bio, let us know directly through the `GitHub Issues <https://github.com/scikit-bio/scikit-bio/issues>`_ page. Your reports are vital for the continuous improvement of scikit-bio.

Wanna contribute? We enthusiastically welcome community contributors! Whether it's adding new features, improving code, or enhancing documentation, your contributions drive scikit-bio and open-source bioinformatics forward. Start your journey by reading the `Contributor's guidelines <https://scikit.bio/contribute.html>`_.


----

Visit the new scikit-bio website: https://scikit.bio to learn more about this project.


Releases
--------

Latest release: `0.6.2 <https://github.com/scikit-bio/scikit-bio/releases/tag/0.6.2>`_ (`documentation <https://scikit.bio/docs/0.6.2/index.html>`_, `changelog <https://github.com/scikit-bio/scikit-bio/blob/main/CHANGELOG.md#version-062>`_). Compatible with Python 3.8 and above.


Installation
------------

Install the latest release of scikit-bio using ``conda``::

    conda install -c conda-forge scikit-bio

Or using ``pip``::

    pip install scikit-bio

Verify the installation::

    python -m skbio.test

See further `instructions on installing <https://scikit.bio/install.html>`_ scikit-bio on various platforms.


Adoption
--------

Some of the projects that we know of that are using scikit-bio are:

- `QIIME 2 <https://qiime2.org/>`_, `Qiita <https://qiita.ucsd.edu/>`_, `Emperor <https://biocore.github.io/emperor/>`_, `tax2tree <https://github.com/biocore/tax2tree>`_, `ghost-tree <https://github.com/JTFouquier/ghost-tree>`_, `Platypus-Conquistador <https://github.com/biocore/Platypus-Conquistador>`_, `An Introduction to Applied Bioinformatics <https://readiab.org>`_.


License
-------

scikit-bio is available under the new BSD license. See `LICENSE.txt <LICENSE.txt>`_ for scikit-bio's license, and the `licenses directory <licenses>`_ for the licenses of third-party software that is (either partially or entirely) distributed with scikit-bio.


Team
----

Our core development team consists of three lead developers: **Dr. Qiyun Zhu** at Arizona State University (ASU) (@qiyunzhu), **Dr. James Morton** at Gutz Analytics (@mortonjt), and **Dr. Daniel McDonald** at the University of California San Diego (UCSD) (@wasade), one software engineer: **Matthew Aton** (@mataton) and one bioinformatician: **Dr. Lars Hunger** (@LarsHunger). **Dr. Rob Knight** at UCSD (@rob-knight) provides guidance on the development and research. **Dr. Greg Caporaso** (@gregcaporaso) at Northern Arizona University (NAU), the former leader of the scikit-bio project, serves as an advisor on the current project.


Credits
-------

We thank the many contributors to scikit-bio. A complete `list of contributors <graphs/contributors>`_ to the scikit-bio codebase is available at GitHub. This however may miss the larger community who contributed by testing the software and providing valuable comments, who we hold equal appreciation to.

Wanna contribute? We enthusiastically welcome community contributors! Whether it's adding new features, improving code, or enhancing documentation, your contributions drive scikit-bio and open-source bioinformatics forward. Start your journey by reading the `Contributor's guidelines <https://scikit.bio/contribute.html>`_.


Funding
-------

The development of scikit-bio is currently supported by the U.S. Department of Energy, Office of Science under award number `DE-SC0024320 <https://genomicscience.energy.gov/compbioawards2023/#Expanding>`_, awarded to Dr. Qiyun Zhu at ASU (lead PI), Dr. James Morton at Gutz Analytics, and Dr. Rob Knight at UCSD.


Citation
--------

If you use scikit-bio for any published research, please see our `Zenodo page <https://zenodo.org/doi/10.5281/zenodo.593387>`_ for how to cite.


Collaboration
-------------

For collaboration inquiries and other formal communications, please reach out to **Dr. Qiyun Zhu** at `qiyun.zhu@asu.edu`. We welcome academic and industrial partnerships to advance our mission.


Branding
--------

The logo of scikit-bio was created by `Alina Prassas <https://cargocollective.com/alinaprassas>`_. Vector and bitmap image files are available at the `logos <logos>`_ directory.


Pre-history
-----------

scikit-bio began from code derived from `PyCogent <https://github.com/pycogent/pycogent>`_ and `QIIME <https://github.com/biocore/qiime>`_, and the contributors and/or copyright holders have agreed to make the code they wrote for PyCogent and/or QIIME available under the BSD license. The contributors to PyCogent and/or QIIME modules that have been ported to scikit-bio are listed below:

- Rob Knight (@rob-knight), Gavin Huttley (@gavinhuttley), Daniel McDonald (@wasade), Micah Hamady, Antonio Gonzalez (@antgonza), Sandra Smit, Greg Caporaso (@gregcaporaso), Jai Ram Rideout (@jairideout), Cathy Lozupone (@clozupone), Mike Robeson (@mikerobeson), Marcin Cieslik, Peter Maxwell, Jeremy Widmann, Zongzhi Liu, Michael Dwan, Logan Knecht (@loganknecht), Andrew Cochran, Jose Carlos Clemente (@cleme), Damien Coy, Levi McCracken, Andrew Butterfield, Will Van Treuren (@wdwvt1), Justin Kuczynski (@justin212k), Jose Antonio Navas Molina (@josenavas), Matthew Wakefield (@genomematt) and Jens Reeder (@jensreeder).


.. |license| image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
   :alt: License
   :target: https://opensource.org/licenses/BSD-3-Clause
.. |build| image:: https://github.com/scikit-bio/scikit-bio/actions/workflows/ci.yml/badge.svg
   :alt: Build Status
   :target: https://github.com/scikit-bio/scikit-bio/actions/workflows/ci.yml
.. |coverage| image:: https://codecov.io/gh/scikit-bio/scikit-bio/graph/badge.svg?token=1qbzC6d2F5 
   :alt: Coverage Status
   :target: https://codecov.io/gh/scikit-bio/scikit-bio
.. |bench| image:: https://img.shields.io/badge/benchmarked%20by-asv-green.svg
   :alt: ASV Benchmarks
   :target: https://s3-us-west-2.amazonaws.com/scikit-bio.org/benchmarks/main/index.html
.. |release| image:: https://img.shields.io/github/v/release/scikit-bio/scikit-bio.svg
   :alt: Release
   :target: https://github.com/scikit-bio/scikit-bio/releases
.. |pypi| image:: https://img.shields.io/pypi/dm/scikit-bio.svg?label=PyPI%20downloads
   :alt: PyPI Downloads
   :target: https://pypi.org/project/scikit-bio/
.. |conda| image:: https://img.shields.io/conda/dn/conda-forge/scikit-bio.svg?label=Conda%20downloads
   :alt: Conda Downloads
   :target: https://anaconda.org/conda-forge/scikit-bio
.. |gitter| image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Gitter
   :target: https://gitter.im/biocore/scikit-bio
