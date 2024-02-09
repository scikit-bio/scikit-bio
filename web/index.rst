.. remove right sidebar

:html_theme.sidebar_secondary.remove:


.. page style and classes

.. raw:: html

   <style>

     .bd-main .bd-content .bd-article-container {
       max-width: 100%;
     }

     details > summary {
       font-weight: bold;
     }

     .centered {
       text-align: center;
     }

     .subtitle {
       text-align: center;
       font-size: 1.2em;
       margin: 20px;
     }

     .heading {
       text-align: center;
       font-size: 1.5em;
       font-weight: 900;
     }

     .no-top-pardding {
       margin-top: 0;
       padding-top: 0;
     }

   </style>


.. hidden page title

.. title:: scikit-bio


.. light/dark logo image

.. image:: _static/img/logo.png
   :class: only-light
   :width: 600 px
   :align: center

.. image:: _static/img/logo_inv.png
   :class: only-dark
   :width: 600 px
   :align: center


.. brief description of the project

.. rst-class:: subtitle

   A community-driven Python library for bioinformatics, providing versatile data structures, algorithms and educational resources.


.. toctree hidden from document but provides header links

.. toctree::
   :hidden:
   :maxdepth: 1

   Install <install>
   Learn <learn>
   Documentation <https://scikit.bio/docs/latest/index.html>
   Contribute <contribute>
   Community <https://github.com/scikit-bio/scikit-bio/discussions>
   Releases <https://github.com/scikit-bio/scikit-bio/blob/master/CHANGELOG.md>
   About <about>


.. grid:: 1 1 1 3

   .. grid-item-card::
      :class-header: centered

      **For Researchers** :octicon:`beaker;1.2em;sd-text-danger`
      ^^^

      Robust, performant and scalable algorithms tailored for the vast landscape of biological data analysis spanning genomics, microbiomics, ecology, evolutionary biology and more. Built to unveil the insights hidden in complex, multi-omic data.
      +++

      .. dropdown example code, implemented using raw html details instead of sphinx-design's dropdown, because the latter has extra edges

      .. raw:: html

         <details>
         <summary>Example</summary>

      .. code-block:: python

         from skbio.tree import TreeNode
         from skbio.diversity import beta_diversity
         from skbio.stats.ordination import pcoa

         data = pd.read_table('data.tsv', index_col=0)
         metadata = pd.read_table('metadata.tsv', index_col=0)
         tree = TreeNode.read('tree.nwk')

         bdiv = beta_diversity(
             'weighted_unifrac', data, ids=data.index, otu_ids=data.columns, tree=tree
         )

         ordi = pcoa(bdiv, number_of_dimensions=3)
         ordi.plot(metadata, column='bodysite')

      .. image:: _static/img/hmp1_pcoa.png
         :alt: PCoA plot

      .. raw:: html

         </details>

   .. grid-item-card::
      :class-header: centered
      
      **For Educators** :octicon:`mortar-board;1.2em;sd-text-info`
      ^^^

      Fundamental bioinformatics algorithms enriched by comprehensive documentation, examples and references, offering a rich resource for classroom and laboratory education (with proven `success <https://readiab.org/>`_). Designed to spark curiosity and foster innovation.
      +++

      .. raw:: html

         <details>
         <summary>Example</summary>

      .. code-block:: python

         from skbio.alignment import global_pairwise_align_protein
         from skbio.sequence.distance import hamming
         from skbio.stats.distance import DistanceMatrix
         from skbio.tree import nj

         def align_dist(seq1, seq2):
             aln = global_pairwise_align_protein(seq1, seq2)[0]
             return hamming(aln[0], aln[1])

         dm = DistanceMatrix.from_iterable(
            seqs, align_dist, keys=ids, validate=False
         )

         tree = nj(dm).root_at_midpoint()
         print(tree.ascii_art())

      ::

                   /-chicken
                  |
         -root----|                    /-rat
                  |          /--------|
                  |         |          \-mouse
                   \--------|
                            |          /-pig
                            |         |
                             \--------|                    /-chimp
                                      |          /--------|
                                       \--------|          \-human
                                                |
                                                 \-monkey


      .. raw:: html

         </details>

   .. grid-item-card::
      :class-header: centered

      **For Developers** :octicon:`git-merge;1.2em;sd-text-warning`
      ^^^

      Industry-standard, production-ready Python codebase featuring a stable, unit-tested API that streamlines development and integration. Licensed under the :repo:`3-Clause BSD <blob/master/LICENSE.txt>`, it provides an expansive platform for both academic research and commercial ventures.
      +++

      .. raw:: html

         <details>
         <summary>Example</summary>

      .. code-block:: python

         def centralize(mat):
             r"""Center data around its geometric average.

             Parameters
             ----------
             mat : array_like, float
                 a matrix of proportions where
                 rows = compositions and
                 columns = components

             Returns
             -------
             numpy.ndarray
                 centered composition matrix

             Examples
             --------
             >>> import numpy as np
             >>> from skbio.stats.composition import centralize
             >>> X = np.array([[.1,.3,.4, .2],[.2,.2,.2,.4]])
             >>> centralize(X)
             array([[ 0.17445763,  0.30216948,  0.34891526,  0.17445763],
                    [ 0.32495488,  0.18761279,  0.16247744,  0.32495488]])

             """
             mat = closure(mat)
             cen = scipy.stats.gmean(mat, axis=0)
             return perturb_inv(mat, cen)

      .. raw:: html

         </details>

----


.. grid:: 1 1 1 2

   .. grid-item::
      :columns: 12 12 5 5
      :padding: 0 0 4 4

      .. rst-class:: heading

         Install

      .. tab-set::
         :class: no-top-pardding

         .. tab-item:: Conda

            .. code-block:: bash

               conda install -c conda-forge scikit-bio

         .. tab-item:: PyPI

            .. code-block:: bash

               pip install scikit-bio

         .. tab-item:: Dev

            .. code-block:: bash

               pip install git+https://github.com/scikit-bio/scikit-bio.git

         .. tab-item:: More

            See detailed :doc:`instructions <install>` on installing scikit-bio on various platforms.

   .. grid-item::
      :columns: 12 12 7 7
      :padding: 0 0 4 4

      .. rst-class:: heading

         News

      .. card-carousel:: 3

         .. card::

            Latest release (:repo:`changelog <blob/master/CHANGELOG.md#version-059>`):

            .. button-link:: https://github.com/scikit-bio/scikit-bio/releases/tag/0.5.9
               :color: success
               :shadow:

               scikit-bio 0.5.9

         .. card::

            New `DOE award <https://genomicscience.energy.gov/compbioawards2023/#Expanding>`_ for scikit-bio development in multi-omics and complex modeling.

         .. card::

            Upcoming scikit-bio `workshop at ISMB 2024 <https://www.iscb.org/ismb2024/programme-schedule/tutorials#ip3>`_, July 11, Montreal, Canada. Welcome to join!

         .. card::

            New website: `scikit.bio <https://scikit.bio>`_ and organization: https://github.com/scikit-bio are online.

----


.. rst-class:: heading

   Feature Highlights

.. grid:: 1 2 3 3
   :padding: 0 0 4 4
   :gutter: 3

   .. grid-item-card::

      :fa:`dna;fa-2x sd-text-success`

      **Biological sequences**: Efficient data structure with a :docs:`flexible grammar <sequence.GrammaredSequence>` for easy manipulation, :docs:`annotation <sequence.GrammaredSequence.has_positional_metadata>`, :docs:`alignment <alignment>`, and conversion into :docs:`motifs <sequence.GrammaredSequence.find_motifs>`, :docs:`k-mers <sequence.GrammaredSequence.iter_kmers>` or tokens for in-depth analysis.

   .. grid-item-card::

      :fa:`tree;fa-2x sd-text-success`

      **Phylogenetic trees**: Scalable :docs:`tree structure <tree.TreeNode>` tailored for evolutionary biology, supporting diverse operations in :docs:`navigation <tree.TreeNode.lca>`, :docs:`manipulation <tree.TreeNode.root_at_midpoint>`, :docs:`comparison <tree.TreeNode.compare_rfd>`, and :docs:`construction <tree.nj>`.

   .. grid-item-card::

      :fa:`chart-column;fa-2x sd-text-success`

      **Community diversity** analysis for ecological studies, with an extensive suite of metrics such as :docs:`UniFrac <diversity.beta.weighted_unifrac>` and :docs:`PD <diversity.alpha.faith_pd>`, optimized to handle large-scale community datasets.

   .. grid-item-card::

      :fa:`compass-drafting;fa-2x sd-text-success`

      **Ordination methods**, such as :docs:`PCoA <stats.ordination.pcoa>`, :docs:`CA <stats.ordination.ca>`, and :docs:`RDA <stats.ordination.rda>`, to uncover patterns underlying high-dimensional data, facilitating insightful visualization.

   .. grid-item-card::

      :fa:`arrow-up-right-dots;fa-2x sd-text-success`

      **Multivariate statistical tests**, such as :docs:`PERMANOVA <stats.distance.permanova>`, :docs:`BIOENV <stats.distance.bioenv>`, and :docs:`Mantel <stats.distance.mantel>`, to decode complex relationships across data matrices and sample properties.

   .. grid-item-card::

      :fa:`chart-pie;fa-2x sd-text-success`

      **Compositional data** processing and analysis, such as :docs:`CLR <stats.composition.clr>` transform and :docs:`ANCOM <stats.composition.ancom>`, built for various omic data types from high-throughput experiments.
