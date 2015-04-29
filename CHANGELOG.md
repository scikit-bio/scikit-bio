# scikit-bio changelog

## Version 0.2.3-dev (changes since 0.2.3 release go here)

### Features
* Added `skbio.stats.composition` for analyzing data made up of proportions
* Added new ``skbio.stats.evolve`` subpackage for evolutionary statistics. Currently contains a single function, ``hommola_cospeciation``, which implements a permutation-based test of correlation between two distance matrices.

### Performance Enhancements
* The speed of quality score decoding has been significantly improved (~2x) when reading `fastq` files.

### Bug fixes
* Changed `BiologicalSequence.distance` to raise an error any time two sequences are passed of different lengths regardless of the `distance_fn` being passed. [(#514)](https://github.com/biocore/scikit-bio/issues/514)
* Fixed issue with ``TreeNode.extend`` where if given the children of another ``TreeNode`` object (``tree.children``), both trees would be left in an incorrect and unpredictable state. ([#889](https://github.com/biocore/scikit-bio/issues/889))

### Deprecated functionality
* Deprecated `skbio.util.flatten`. This function will be removed in scikit-bio 0.3.1. Please use standard python library functionality
described here [Making a flat list out of lists of lists](http://stackoverflow.com/a/952952/3639023), [Flattening a shallow list](http://stackoverflow.com/a/406199/3639023) ([#833](https://github.com/biocore/scikit-bio/issues/833))

### Backward-incompatible changes
* Removed the following deprecated functionality:
    - `skbio.parse` subpackage, including `SequenceIterator`, `FastaIterator`, `FastqIterator`, `load`, `parse_fasta`, `parse_fastq`, `parse_qual`, `write_clustal`, `parse_clustal`, and `FastqParseError`; please use `skbio.io` instead.
    - `skbio.format` subpackage, including `fasta_from_sequence`, `fasta_from_alignment`, and `format_fastq_record`; please use `skbio.io` instead.
    - `skbio.alignment.SequenceCollection.int_map`; please use `SequenceCollection.update_ids` instead.
    - `skbio.alignment.SequenceCollection` methods `to_fasta` and `toFasta`; please use `SequenceCollection.write` instead.
    - `constructor` parameter in `skbio.alignment.Alignment.majority_consensus`; please convert returned biological sequence object manually as desired (e.g., `str(seq)`).
    - `skbio.alignment.Alignment.to_phylip`; please use `Alignment.write` instead.
    - `skbio.sequence.BiologicalSequence.to_fasta`; please use `BiologicalSequence.write` instead.
    - `skbio.tree.TreeNode` methods `from_newick`, `from_file`, and `to_newick`; please use `TreeNode.read` and `TreeNode.write` instead.
    - `skbio.stats.distance.DissimilarityMatrix` methods `from_file` and `to_file`; please use `DissimilarityMatrix.read` and `DissimilarityMatrix.write` instead.
    - `skbio.stats.ordination.OrdinationResults` methods `from_file` and `to_file`; please use `OrdinationResults.read` and `OrdinationResults.write` instead.
    - `skbio.stats.p_value_to_str`; there is no replacement.
    - `skbio.stats.subsample`; please use `skbio.stats.subsample_counts` instead.
    - `skbio.stats.distance.ANOSIM`; please use `skbio.stats.distance.anosim` instead.
    - `skbio.stats.distance.PERMANOVA`; please use `skbio.stats.distance.permanova` instead.
    - `skbio.stats.distance.CategoricalStatsResults`; there is no replacement, please use `skbio.stats.distance.anosim` or `skbio.stats.distance.permanova`, which will return a `pandas.Series` object.
* `skbio.alignment.Alignment.majority_consensus` now returns `BiologicalSequence('')` if the alignment is empty. Previously, `''` was returned.

### Miscellaneous
* Changed testing procedures
    - Developers should now use `make test`
    - Users can use `python -m skbio.test`
    - Added `skbio.util._testing.TestRunner` (available through `skbio.util.TestRunner`). Used to provide a `test` method for each module init file. This class represents a unified testing path which wraps all `skbio` testing functionality.
    - Autodetect Python version and disable doctests for Python 3.
* `numpy` is no longer required to be installed before installing scikit-bio!
* Upgraded checklist.py to check source files non-conforming to [new header style](http://scikit-bio.org/docs/latest/development/new_module.html). ([#855](https://github.com/biocore/scikit-bio/issues/855))

## Version 0.2.3 (2015-02-13)

### Features
* Modified ``skbio.stats.distance.pwmantel`` to accept a list of filepaths. This is useful as it allows for a smaller amount of memory consumption as it only loads two matrices at a time as opposed to requiring that all distance matrices are loaded into memory.
* Added ``skbio.util.find_duplicates`` for finding duplicate elements in an iterable.

### Bug fixes
* Fixed floating point precision bugs in ``Alignment.position_frequencies``, ``Alignment.position_entropies``, ``Alignment.omit_gap_positions``, ``Alignment.omit_gap_sequences``, ``BiologicalSequence.k_word_frequencies``, and ``SequenceCollection.k_word_frequencies`` ([#801](https://github.com/biocore/scikit-bio/issues/801)).

### Backward-incompatible changes
* Removed ``feature_types`` attribute from ``BiologicalSequence`` and all subclasses ([#797](https://github.com/biocore/scikit-bio/pull/797)).
* Removed ``find_features`` method from ``BiologicalSequence`` and ``ProteinSequence`` ([#797](https://github.com/biocore/scikit-bio/pull/797)).
* ``BiologicalSequence.k_word_frequencies`` now returns a ``collections.defaultdict`` of type ``float`` instead of type ``int``. This only affects the "default" case, when a key isn't present in the dictionary. Previous behavior would return ``0`` as an ``int``, while the new behavior is to return ``0.0`` as a ``float``. This change also affects the ``defaultdict``s that are returned by ``SequenceCollection.k_word_frequencies``.

### Miscellaneous
* ``DissimilarityMatrix`` and ``DistanceMatrix`` now report duplicate IDs in the ``DissimilarityMatrixError`` message that can be raised during validation.

## Version 0.2.2 (2014-12-04)

### Features
* Added ``plot`` method to ``skbio.stats.distance.DissimilarityMatrix`` for creating basic heatmaps of a dissimilarity/distance matrix (see [#684](https://github.com/biocore/scikit-bio/issues/684)). Also added  ``_repr_png_`` and ``_repr_svg_`` methods for automatic display in the IPython Notebook, with ``png`` and ``svg`` properties for direct access.
* Added `__str__` method to `skbio.stats.ordination.OrdinationResults`.
* Added ``skbio.stats.distance.anosim`` and ``skbio.stats.distance.permanova`` functions, which replace the ``skbio.stats.distance.ANOSIM`` and ``skbio.stats.distance.PERMANOVA`` classes. These new functions provide simpler procedural interfaces to running these statistical methods. They also provide more convenient access to results by returning a ``pandas.Series`` instead of a ``CategoricalStatsResults`` object. These functions have more extensive documentation than their previous versions. If significance tests are suppressed, p-values are returned as ``np.nan`` instead of ``None`` for consistency with other statistical methods in scikit-bio. [#754](https://github.com/biocore/scikit-bio/issues/754)
* Added `skbio.stats.power` for performing empirical power analysis. The module uses existing datasets and iteratively draws samples to estimate the number of samples needed to see a significant difference for a given critical value.
* Added `skbio.stats.isubsample` for subsampling from an unknown number of values. This method supports subsampling from multiple partitions and does not require that all items be stored in memory, requiring approximately `O(N*M)`` space where `N` is the number of partitions and `M` is the maximum subsample size.
* Added ``skbio.stats.subsample_counts``, which replaces ``skbio.stats.subsample``. See deprecation section below for more details ([#770](https://github.com/biocore/scikit-bio/issues/770)).

### Bug fixes
* Fixed issue where SSW wouldn't compile on i686 architectures ([#409](https://github.com/biocore/scikit-bio/issues/409)).

### Deprecated functionality
* Deprecated ``skbio.stats.p_value_to_str``. This function will be removed in scikit-bio 0.3.0. Permutation-based p-values in scikit-bio are calculated as ``(num_extreme + 1) / (num_permutations + 1)``, so it is impossible to obtain a p-value of zero. This function historically existed for correcting the number of digits displayed when obtaining a p-value of zero. Since this is no longer possible, this functionality will be removed.
* Deprecated ``skbio.stats.distance.ANOSIM`` and ``skbio.stats.distance.PERMANOVA`` in favor of ``skbio.stats.distance.anosim`` and ``skbio.stats.distance.permanova``, respectively.
* Deprecated ``skbio.stats.distance.CategoricalStatsResults`` in favor of using ``pandas.Series`` to store statistical method results. ``anosim`` and ``permanova`` return ``pandas.Series`` instead of ``CategoricalStatsResults``.
* Deprecated ``skbio.stats.subsample`` in favor of ``skbio.stats.subsample_counts``, which provides an identical interface; only the function name has changed. ``skbio.stats.subsample`` will be removed in scikit-bio 0.3.0.

### Backward-incompatible changes
* Deprecation warnings are now raised using ``DeprecationWarning`` instead of ``UserWarning`` ([#774](https://github.com/biocore/scikit-bio/issues/774)).

### Miscellaneous
* The ``pandas.DataFrame`` returned by ``skbio.stats.distance.pwmantel`` now stores p-values as floats and does not convert them to strings with a specific number of digits. p-values that were previously stored as "N/A" are now stored as ``np.nan`` for consistency with other statistical methods in scikit-bio. See note in "Deprecated functionality" above regarding ``p_value_to_str`` for details.
* scikit-bio now supports versions of IPython < 2.0.0 ([#767](https://github.com/biocore/scikit-bio/issues/767)).

## Version 0.2.1 (2014-10-27)

This is an alpha release of scikit-bio. At this stage, major backwards-incompatible API changes can and will happen. Unified I/O with the scikit-bio I/O registry was the focus of this release.

### Features
* Added ``strict`` and ``lookup`` optional parameters to ``skbio.stats.distance.mantel`` for handling reordering and matching of IDs when provided ``DistanceMatrix`` instances as input (these parameters were previously only available in ``skbio.stats.distance.pwmantel``).
* ``skbio.stats.distance.pwmantel`` now accepts an iterable of ``array_like`` objects. Previously, only ``DistanceMatrix`` instances were allowed.
* Added ``plot`` method to ``skbio.stats.ordination.OrdinationResults`` for creating basic 3-D matplotlib scatterplots of ordination results, optionally colored by metadata in a ``pandas.DataFrame`` (see [#518](https://github.com/biocore/scikit-bio/issues/518)). Also added  ``_repr_png_`` and ``_repr_svg_`` methods for automatic display in the IPython Notebook, with ``png`` and ``svg`` properties for direct access.
* Added ``skbio.stats.ordination.assert_ordination_results_equal`` for comparing ``OrdinationResults`` objects for equality in unit tests.
* ``BiologicalSequence`` (and its subclasses) now optionally store Phred quality scores. A biological sequence's quality scores are stored as a 1-D ``numpy.ndarray`` of nonnegative integers that is the same length as the biological sequence. Quality scores can be provided upon object instantiation via the keyword argument ``quality``, and can be retrieved via the ``BiologicalSequence.quality`` property. ``BiologicalSequence.has_quality`` is also provided for determining whether a biological sequence has quality scores or not. See [#616](https://github.com/biocore/scikit-bio/issues/616) for more details.
* Added ``BiologicalSequence.sequence`` property for retrieving the underlying string representing the sequence characters. This was previously (and still is) accessible via ``BiologicalSequence.__str__``. It is provided via a property for convenience and explicitness.
* Added ``BiologicalSequence.equals`` for full control over equality testing of biological sequences. By default, biological sequences must have the same type, underlying sequence of characters, identifier, description, and quality scores to compare equal. These properties can be ignored via the keyword argument ``ignore``. The behavior of ``BiologicalSequence.__eq__``/``__ne__`` remains unchanged (only type and underlying sequence of characters are compared).
* Added ``BiologicalSequence.copy`` for creating a copy of a biological sequence, optionally with one or more attributes updated.
* ``BiologicalSequence.__getitem__`` now supports specifying a sequence of indices to take from the biological sequence.
* Methods to read and write taxonomies are now available under ``skbio.tree.TreeNode.from_taxonomy`` and ``skbio.tree.TreeNode.to_taxonomy`` respectively.
* Added ``SequenceCollection.update_ids``, which provides a flexible way of updating sequence IDs on a ``SequenceCollection`` or ``Alignment`` (note that a new object is returned, since instances of these classes are immutable). Deprecated ``SequenceCollection.int_map`` in favor of this new method; it will be removed in scikit-bio 0.3.0.
* Added ``skbio.util.cardinal_to_ordinal`` for converting a cardinal number to ordinal string (e.g., useful for error messages).
* New I/O Registry: supports multiple file formats, automatic file format detection when reading, unified procedural ``skbio.io.read`` and ``skbio.io.write`` in addition to OOP interfaces (``read/write`` methods) on the below objects. See ``skbio.io`` for more details.
    - Added "clustal" format support:
        * Has sniffer
        * Readers: ``Alignment``
        * Writers: ``Alignment``
    - Added "lsmat" format support:
        * Has sniffer
        * Readers: ``DissimilarityMatrix``, ``DistanceMatrix``
        * Writers: ``DissimilarityMatrix``, ``DistanceMatrix``
    - Added "ordination" format support:
        * Has sniffer
        * Readers: ``OrdinationResults``
        * Writers: ``OrdinationResults``
    - Added "newick" format support:
        * Has sniffer
        * Readers: ``TreeNode``
        * Writers: ``TreeNode``
    - Added "phylip" format support:
        * No sniffer
        * Readers: None
        * Writers: ``Alignment``
    - Added "qseq" format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: None
    - Added "fasta"/QUAL format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``Alignment``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: same as readers
    - Added "fastq" format support:
        * Has sniffer
        * Readers: generator of ``BiologicalSequence`` or its subclasses, ``SequenceCollection``, ``Alignment``, ``BiologicalSequence``, ``NucleotideSequence``, ``DNASequence``, ``RNASequence``, ``ProteinSequence``
        * Writers: same as readers

### Bug fixes

* Removed ``constructor`` parameter from ``Alignment.k_word_frequencies``, ``BiologicalSequence.k_words``, ``BiologicalSequence.k_word_counts``, and ``BiologicalSequence.k_word_frequencies`` as it had no effect (it was never hooked up in the underlying code). ``BiologicalSequence.k_words`` now returns a generator of ``BiologicalSequence`` objects instead of strings.
* Modified the ``Alignment`` constructor to verify that all sequences have the same length, if not, raise an ``AlignmentError`` exception.  Updated the method ``Alignment.subalignment`` to calculate the indices only once now that identical sequence length is guaranteed.

### Deprecated functionality
* Deprecated ``constructor`` parameter in ``Alignment.majority_consensus`` in favor of having users call ``str`` on the returned ``BiologicalSequence``. This parameter will be removed in scikit-bio 0.3.0.

* Existing I/O functionality deprecated in favor of I/O registry, old functionality will be removed in scikit-bio 0.3.0. All functionality can be found at ``skbio.io.read``, ``skbio.io.write``, and the methods listed below:
    * Deprecated the following "clustal" readers/writers:
        - ``write_clustal`` -> ``Alignment.write``
        - ``parse_clustal`` -> ``Alignment.read``

    * Deprecated the following distance matrix format ("lsmat") readers/writers:
        - ``DissimilarityMatrix.from_file`` -> ``DissimilarityMatrix.read``
        - ``DissimilarityMatrix.to_file`` -> ``DissimilarityMatrix.write``
        - ``DistanceMatrix.from_file`` -> ``DistanceMatrix.read``
        - ``DistanceMatrix.to_file`` -> ``DistanceMatrix.write``

    * Deprecated the following ordination format ("ordination") readers/writers:
        - ``OrdinationResults.from_file`` -> ``OrdinationResults.read``
        - ``OrdinationResults.to_file`` -> ``OrdinationResults.write``

    * Deprecated the following "newick" readers/writers:
        - ``TreeNode.from_file`` -> ``TreeNode.read``
        - ``TreeNode.from_newick`` -> ``TreeNode.read``
        - ``TreeNode.to_newick`` -> ``TreeNode.write``

    * Deprecated the following "phylip" writers:
        - ``Alignment.to_phylip`` -> ``Alignment.write``

    * Deprecated the following "fasta"/QUAL readers/writers:
        - ``SequenceCollection.from_fasta_records`` -> ``SequenceCollection.read``
        - ``SequenceCollection.to_fasta`` -> ``SequenceCollection.write``
        - ``fasta_from_sequences`` -> ``skbio.io.write(obj, into=<file>, format='fasta')``
        - ``fasta_from_alignment`` -> ``Alignment.write``
        - ``parse_fasta`` -> ``skbio.io.read(<fasta>, format='fasta')``
        - ``parse_qual`` -> ``skbio.io.read(<fasta>, format='fasta', qual=<file>)``
        - ``BiologicalSequence.to_fasta`` -> ``BiologicalSequence.write``

    * Deprecated the following "fastq" readers/writers:
        - ``parse_fastq`` -> ``skbio.io.read(<fastq>, format='fastq')``
        - ``format_fastq_record`` -> ``skbio.io.write(<fastq>, format='fastq')``

### Backward-incompatible changes

* ``skbio.stats.distance.mantel`` now returns a 3-element tuple containing correlation coefficient, p-value, and the number of matching rows/cols in the distance matrices (``n``). The return value was previously a 2-element tuple containing only the correlation coefficient and p-value.
* ``skbio.stats.distance.mantel`` reorders input ``DistanceMatrix`` instances based on matching IDs (see optional parameters ``strict`` and ``lookup`` for controlling this behavior). In the past, ``DistanceMatrix`` instances were treated the same as ``array_like`` input and no reordering took place, regardless of ID (mis)matches. ``array_like`` input behavior remains the same.
* If mismatched types are provided to ``skbio.stats.distance.mantel`` (e.g., a ``DistanceMatrix`` and ``array_like``), a ``TypeError`` will be raised.

### Miscellaneous

* Added git timestamp checking to checklist.py, ensuring that when changes are made to Cython (.pyx) files, their corresponding generated C files are also updated.
* Fixed performance bug when instantiating ``BiologicalSequence`` objects. The previous runtime scaled linearly with sequence length; it is now constant time when the sequence is already a string. See [#623](https://github.com/biocore/scikit-bio/issues/623) for details.
* IPython and six are now required dependencies.

## Version 0.2.0 (2014-08-07)

This is an initial alpha release of scikit-bio. At this stage, major backwards-incompatible API changes can and will happen. Many backwards-incompatible API changes were made since the previous release.

### Features

* Added ability to compute distances between sequences in a ``SequenceCollection`` object ([#509](https://github.com/biocore/scikit-bio/issues/509)), and expanded ``Alignment.distance`` to allow the user to pass a function for computing distances (the default distance metric is still ``scipy.spatial.distance.hamming``) ([#194](https://github.com/biocore/scikit-bio/issues/194)).
* Added functionality to not penalize terminal gaps in global alignment. This functionality results in more biologically relevant global alignments (see [#537](https://github.com/biocore/scikit-bio/issues/537) for discussion of the issue) and is now the default behavior for global alignment.
* The python global aligners (``global_pairwise_align``, ``global_pairwise_align_nucleotide``, and ``global_pairwise_align_protein``) now support aligning pairs of sequences, pairs of alignments, and a sequence and an alignment (see [#550](https://github.com/biocore/scikit-bio/issues/550)). This functionality supports progressive multiple sequence alignment, among other things such as adding a sequence to an existing alignment.
* Added ``StockholmAlignment.to_file`` for writing Stockholm-formatted files.
* Added ``strict=True`` optional parameter to ``DissimilarityMatrix.filter``.
* Added ``TreeNode.find_all`` for finding all tree nodes that match a given name.


### Bug fixes

* Fixed bug that resulted in a ``ValueError`` from ``local_align_pairwise_nucleotide`` (see [#504](https://github.com/biocore/scikit-bio/issues/504)) under many circumstances. This would not generate incorrect results, but would cause the code to fail.

### Backward-incompatible changes

* Removed ``skbio.math``, leaving ``stats`` and ``diversity`` to become top level packages. For example, instead of ``from skbio.math.stats.ordination import PCoA`` you would now import ``from skbio.stats.ordination import PCoA``.
* The module ``skbio.math.gradient`` as well as the contents of ``skbio.math.subsample`` and ``skbio.math.stats.misc`` are now found in ``skbio.stats``. As an example, to import subsample: ``from skbio.stats import subsample``; to import everything from gradient: ``from skbio.stats.gradient import *``.
* The contents of ``skbio.math.stats.ordination.utils`` are now in ``skbio.stats.ordination``.
* Removed ``skbio.app`` subpackage (i.e., the *application controller framework*) as this code has been ported to the standalone [burrito](https://github.com/biocore/burrito) Python package. This code was not specific to bioinformatics and is useful for wrapping command-line applications in general.
* Removed ``skbio.core``, leaving ``alignment``, ``genetic_code``, ``sequence``, ``tree``, and ``workflow`` to become top level packages. For example, instead of ``from skbio.core.sequence import DNA`` you would now import ``from skbio.sequence import DNA``.
* Removed ``skbio.util.exception`` and ``skbio.util.warning`` (see [#577](https://github.com/biocore/scikit-bio/issues/577) for the reasoning behind this change). The exceptions/warnings were moved to the following locations:
 - ``FileFormatError``, ``RecordError``, ``FieldError``, and ``EfficiencyWarning`` have been moved to ``skbio.util``
 - ``BiologicalSequenceError`` has been moved to ``skbio.sequence``
 - ``SequenceCollectionError`` and ``StockholmParseError`` have been moved to ``skbio.alignment``
 - ``DissimilarityMatrixError``, ``DistanceMatrixError``, ``DissimilarityMatrixFormatError``, and ``MissingIDError`` have been moved to ``skbio.stats.distance``
 - ``TreeError``, ``NoLengthError``, ``DuplicateNodeError``, ``MissingNodeError``, and ``NoParentError`` have been moved to ``skbio.tree``
 - ``FastqParseError`` has been moved to ``skbio.parse.sequences``
 - ``GeneticCodeError``, ``GeneticCodeInitError``, and ``InvalidCodonError`` have been moved to ``skbio.genetic_code``
* The contents of ``skbio.genetic_code`` formerly ``skbio.core.genetic_code`` are now in ``skbio.sequence``. The ``GeneticCodes`` dictionary is now a function ``genetic_code``. The functionality is the same, except that because this is now a function rather than a dict, retrieving a genetic code is done using a function call rather than a lookup (so, for example, ``GeneticCodes[2]`` becomes ``genetic_code(2)``.
* Many submodules have been made private with the intention of simplifying imports for users. See [#562](https://github.com/biocore/scikit-bio/issues/562) for discussion of this change. The following list contains the previous module name and where imports from that module should now come from.
 - ``skbio.alignment.ssw`` to ``skbio.alignment``
 - ``skbio.alignment.alignment`` to ``skbio.alignment``
 - ``skbio.alignment.pairwise`` to ``skbio.alignment``
 - ``skbio.diversity.alpha.base`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.alpha.gini`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.alpha.lladser`` to ``skbio.diversity.alpha``
 - ``skbio.diversity.beta.base`` to ``skbio.diversity.beta``
 - ``skbio.draw.distributions`` to ``skbio.draw``
 - ``skbio.stats.distance.anosim`` to ``skbio.stats.distance``
 - ``skbio.stats.distance.base`` to ``skbio.stats.distance``
 - ``skbio.stats.distance.permanova`` to ``skbio.stats.distance``
 - ``skbio.distance`` to ``skbio.stats.distance``
 - ``skbio.stats.ordination.base`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.canonical_correspondence_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.correspondence_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.principal_coordinate_analysis`` to ``skbio.stats.ordination``
 - ``skbio.stats.ordination.redundancy_analysis`` to ``skbio.stats.ordination``
 - ``skbio.tree.tree`` to ``skbio.tree``
 - ``skbio.tree.trie`` to ``skbio.tree``
 - ``skbio.util.misc`` to ``skbio.util``
 - ``skbio.util.testing`` to ``skbio.util``
 - ``skbio.util.exception`` to ``skbio.util``
 - ``skbio.util.warning`` to ``skbio.util``
* Moved ``skbio.distance`` contents into ``skbio.stats.distance``.

### Miscellaneous

* Relaxed requirement in ``BiologicalSequence.distance`` that sequences being compared are of equal length. This is relevant for Hamming distance, so the check is still performed in that case, but other distance metrics may not have that requirement. See [#504](https://github.com/biocore/scikit-bio/issues/507)).
* Renamed ``powertrip.py`` repo-checking script to ``checklist.py`` for clarity.
* ``checklist.py`` now ensures that all unit tests import from a minimally deep API. For example, it will produce an error if ``skbio.core.distance.DistanceMatrix`` is used over ``skbio.DistanceMatrix``.
* Extra dimension is no longer calculated in ``skbio.stats.spatial.procrustes``.
* Expanded documentation in various subpackages.
* Added new scikit-bio logo. Thanks [Alina Prassas](http://cargocollective.com/alinaprassas)!

## Version 0.1.4 (2014-06-25)

This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

### Features

* Added Python implementations of Smith-Waterman and Needleman-Wunsch alignment as ``skbio.core.alignment.pairwise.local_pairwise_align`` and ``skbio.core.alignment.pairwise.global_pairwise_align``. These are much slower than native C implementations (e.g., ``skbio.core.alignment.local_pairwise_align_ssw``) and as a result raise an ``EfficencyWarning`` when called, but are included as they serve as useful educational examples as they’re simple to experiment with.
* Added ``skbio.core.diversity.beta.pw_distances`` and ``skbio.core.diversity.beta.pw_distances_from_table``. These provide convenient access to the ``scipy.spatial.distance.pdist`` *beta diversity* metrics from within scikit-bio. The ``skbio.core.diversity.beta.pw_distances_from_table`` function will only be available temporarily, until the ``biom.table.Table`` object is merged into scikit-bio (see [#489](https://github.com/biocore/scikit-bio/issues/489)), at which point ``skbio.core.diversity.beta.pw_distances`` will be updated to use that.
* Added ``skbio.core.alignment.StockholmAlignment``, which provides support for parsing [Stockholm-formatted alignment files](http://sonnhammer.sbc.su.se/Stockholm.html) and working with those alignments in the context RNA secondary structural information.
* Added ``skbio.core.tree.majority_rule`` function for computing consensus trees from a list of trees.

### Backward-incompatible changes

* Function ``skbio.core.alignment.align_striped_smith_waterman`` renamed to ``local_pairwise_align_ssw`` and now returns an ``Alignment`` object instead of an ``AlignmentStructure``
* The following keyword-arguments for ``StripedSmithWaterman`` and ``local_pairwise_align_ssw`` have been renamed:
    * ``gap_open`` -> ``gap_open_penalty``
    * ``gap_extend`` -> ``gap_extend_penalty``
    * ``match`` -> ``match_score``
    * ``mismatch`` -> ``mismatch_score``
* Removed ``skbio.util.sort`` module in favor of [natsort](https://pypi.python.org/pypi/natsort) package.

### Miscellaneous

* Added powertrip.py script to perform basic sanity-checking of the repo based on recurring issues that weren't being caught until release time; added to Travis build.
* Added RELEASE.md with release instructions.
* Added intersphinx mappings to docs so that "See Also" references to numpy, scipy, matplotlib, and pandas are hyperlinks.
* The following classes are no longer ``namedtuple`` subclasses (see [#359](https://github.com/biocore/scikit-bio/issues/359) for the rationale):
    * ``skbio.math.stats.ordination.OrdinationResults``
    * ``skbio.math.gradient.GroupResults``
    * ``skbio.math.gradient.CategoryResults``
    * ``skbio.math.gradient.GradientANOVAResults``
* Added coding guidelines draft.
* Added new alpha diversity formulas to the ``skbio.math.diversity.alpha`` documentation.

## Version 0.1.3 (2014-06-12)

This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

### Features

* Added ``enforce_qual_range`` parameter to ``parse_fastq`` (on by default, maintaining backward compatibility). This allows disabling of the quality score range-checking.
* Added ``skbio.core.tree.nj``, which applies neighbor-joining for phylogenetic reconstruction.
* Added ``bioenv``, ``mantel``, and ``pwmantel`` distance-based statistics to ``skbio.math.stats.distance`` subpackage.
* Added ``skbio.math.stats.misc`` module for miscellaneous stats utility functions.
* IDs are now optional when constructing a ``DissimilarityMatrix`` or ``DistanceMatrix`` (monotonically-increasing integers cast as strings are automatically used).
* Added ``DistanceMatrix.permute`` method for randomly permuting rows and columns of a distance matrix.
* Added the following methods to ``DissimilarityMatrix``: ``filter``, ``index``, and ``__contains__`` for ID-based filtering, index lookup, and membership testing, respectively.
* Added ``ignore_comment`` parameter to ``parse_fasta`` (off by default, maintaining backward compatibility). This handles stripping the comment field from the header line (i.e., all characters beginning with the first space) before returning the label.
* Added imports of ``BiologicalSequence``, ``NucleotideSequence``, ``DNA``, ``DNASequence``, ``RNA``, ``RNASequence``, ``Protein``, ``ProteinSequence``, ``DistanceMatrix``, ``align_striped_smith_waterman``, `` SequenceCollection``, ``Alignment``, ``TreeNode``, ``nj``, ``parse_fasta``, ``parse_fastq``, ``parse_qual``, ``FastaIterator``, ``FastqIterator``, ``SequenceIterator`` in ``skbio/__init__.py`` for convenient importing. For example, it's now possible to ``from skbio import Alignment``, rather than ``from skbio.core.alignment import Alignment``.

### Bug fixes

* Fixed a couple of unit tests that could fail stochastically.
* Added missing ``__init__.py`` files to a couple of test directories so that these tests won't be skipped.
* ``parse_fastq`` now raises an error on dangling records.
* Fixed several warnings that were raised while running the test suite with Python 3.4.

### Backward-incompatible changes

* Functionality imported from ``skbio.core.ssw`` must now be imported from ``skbio.core.alignment`` instead.

### Miscellaneous

* Code is now flake8-compliant; added flake8 checking to Travis build.
* Various additions and improvements to documentation (API, installation instructions, developer instructions, etc.).
* ``__future__`` imports are now standardized across the codebase.
* New website front page and styling changes throughout. Moved docs site to its own versioned subdirectories.
* Reorganized alignment data structures and algorithms (e.g., SSW code, ``Alignment`` class, etc.) into an ``skbio.core.alignment`` subpackage.

## Version 0.1.1 (2014-05-16)

Fixes to setup.py. This is a pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.

## Version 0.1.0 (2014-05-15)

Initial pre-alpha release. At this stage, major backwards-incompatible API changes can and will happen.
