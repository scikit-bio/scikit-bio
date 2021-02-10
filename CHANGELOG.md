# scikit-bio changelog

## Version 0.5.6-dev

### Features

* Added support for `np.float32` with `DissimilarityMatrix` objects.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

* Avoid an implicit data copy on construction of `DissimilarityMatrix` objects.
* Use a memory-optimized version of permute in DistanceMatrix, see [PR #1756](https://github.com/biocore/scikit-bio/pull/1756).
* Refactor pearson and spearman skbio.stats.distance.mantel implementations to drastically improve memory locality. Also cache intermediate results that are invariant across permutations, see [PR #1756](https://github.com/biocore/scikit-bio/pull/1756).

### Bug fixes

* Fix windows and 32bit incompatibility in `unweighted_unifrac`.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* Specify build dependencies in pyproject.toml. This allows the package to be installed without having to first manually install numpy.

* Update hdmedians package to a version which doesn't require an initial manual numpy install.

## Version 0.5.6

### Features

* Added option to return a capture group compiled regex pattern to any class inheriting ``GrammaredSequence`` through the ``to_regex`` method. ([#1431](https://github.com/biocore/scikit-bio/issues/1431))

* Added `Dissimilarity.within` and `.between` to obtain the respective distances and express them as a `DataFrame`. ([#1662](https://github.com/biocore/scikit-bio/pull/1662))

* Added Kendall Tau as possible correlation method in the `skbio.stats.distance.mantel` function ([#1675](https://github.com/biocore/scikit-bio/issues/1675)).

* Added support for IUPAC amino acid codes U (selenocysteine), O (pyrrolysine), and J (leucine or isoleucine). ([#1576](https://github.com/biocore/scikit-bio/issues/1576)

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

* Changed `skbio.tree.TreeNode.support` from a method to a property.
* Added `assign_supports` method to `skbio.tree.TreeNode` to extract branch support values from node labels.
* Modified the way a node's label is printed: `support:name` if both exist, or `support` or `name` if either exists.

### Performance enhancements

### Bug fixes
* Require `Sphinx <= 3.0`. Newer Sphinx versions caused build errors. [#1719](https://github.com/biocore/scikit-bio/pull/1719)

* * `skbio.stats.ordination` tests have been relaxed. ([#1713](https://github.com/biocore/scikit-bio/issues/1713))

* Fixes build errors for newer versions of NumPy, Pandas, and SciPy.

* Corrected a criticial bug in `skbio.alignment.StripedSmithWaterman`/`skbio.alignment.local_pairwise_align_ssw` which would cause the formatting of the aligned sequences to misplace gap characters by the number of gap characters present in the opposing aligned sequence up to that point. This was caused by a faulty implementation of CIGAR string parsing, see [#1679](https://github.com/biocore/scikit-bio/pull/1679) for full details.

* Fixes build errors for newer versions of NumPy, Pandas, and SciPy.

* Corrected a criticial bug in `skbio.alignment.StripedSmithWaterman`/`skbio.alignment.local_pairwise_align_ssw` which would cause the formatting of the aligned sequences to misplace gap characters by the number of gap characters present in the opposing aligned sequence up to that point. This was caused by a faulty implementation of CIGAR string parsing, see [#1679](https://github.com/biocore/scikit-bio/pull/1679) for full details.

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

* `skbio.diversity.beta_diversity` now accepts a pandas DataFrame as input.

* Avoid pandas 1.0.0 import warning ([#1688](https://github.com/biocore/scikit-bio/issues/1688))

* Added support for Python 3.8 and dropped support for Python 3.5.

* This version now depends on `scipy >= 1.3` and `pandas >= 1.0`.

## Version 0.5.5 (2018-12-10)

### Features

* `skbio.stats.composition` now has methods to compute additive log-ratio transformation and inverse additive log-ratio transformation (`alr`, `alr_inv`) as well as a method to build a basis from a sequential binary partition (`sbp_basis`).

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous
* Python 3.6 and 3.7 compatibility is now supported

* A pytest runner is shipped with every installation ([#1633](https://github.com/biocore/scikit-bio/pull/1633))

* The nosetest framework has been replaced in favor of pytest ([#1624](https://github.com/biocore/scikit-bio/pull/1624))

* The numpy docs are deprecated in favor of [Napoleon](http://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html) ([#1629](https://github.com/biocore/scikit-bio/pull/1629))

* This version is now compatible with NumPy >= 1.9.2 and Pandas >= 0.23. ([#1627](https://github.com/biocore/scikit-bio/pull/1627))

## Version 0.5.4 (2018-08-23)

### Features

* Added `FSVD`, an alternative fast heuristic method to perform Principal Coordinates Analysis, to `skbio.stats.ordination.pcoa`.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

* Added optimized utility methods `f_matrix_inplace` and `e_matrix_inplace` which perform `f_matrix` and `e_matrix` computations in-place and are used by the new `center_distance_matrix` method in `skbio.stats.ordination`.

### Bug fixes

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous

## Version 0.5.3 (2018-08-07)

### Features
* Added `unpack` and `unpack_by_func` methods to `skbio.tree.TreeNode` to unpack one or multiple internal nodes. The `unpack` operation removes an internal node and regrafts its children to its parent while retaining the overall length. ([#1572](https://github.com/biocore/scikit-bio/pull/1572))
* Added `support` to `skbio.tree.TreeNode` to return the support value of a node.
* Added `permdisp` to `skbio.stats.distance` to test for the homogeniety of groups. ([#1228](https://github.com/biocore/scikit-bio/issues/1228)).

* Added `pcoa_biplot` to `skbio.stats.ordination` to project descriptors into a PCoA plot.

* Fixed pandas to 0.22.0 due to this: https://github.com/pandas-dev/pandas/issues/20527

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements

### Bug fixes

* Relaxing type checking in diversity calculations.  ([#1583](https://github.com/biocore/scikit-bio/issues/1583)).

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous


## Version 0.5.2 (2018-04-18)

### Features
* Added ``skbio.io.format.embl`` for reading and writing EMBL files for ``DNA``, ``RNA`` and ``Sequence`` classes.

* Removing ValueError check in `skbio.stats._subsample.subsample_counts` when `replace=True` and `n` is greater than the number of items in counts.  [#1527](https://github.com/biocore/scikit-bio/pull/1527)

* Added ``skbio.io.format.gff3`` for reading and writing GFF3 files for ``DNA``, ``Sequence``, and ``IntervalMetadata`` classes. ([#1450](https://github.com/biocore/scikit-bio/pull/1450))

* `skbio.metadata.IntervalMetadata` constructor has a new keyword argument, `copy_from`, for creating an `IntervalMetadata` object from an existing `IntervalMetadata` object with specified `upper_bound`.

* `skbio.metadata.IntervalMetadata` constructor allows `None` as a valid value for `upper_bound`. An `upper_bound` of `None` means that the `IntervalMetadata` object has no upper bound.

* `skbio.metadata.IntervalMetadata.drop` has a new boolean parameter `negate` to indicate whether to drop or keep the specified `Interval` objects.

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]

### Performance enhancements
* `skbio.tree.nj` wall-clock runtime was decreased by 99% for a 500x500 distance matrix and 93% for a 100x100 distance matrix. ([#1512](https://github.com/biocore/scikit-bio/pull/1512), [#1513](https://github.com/biocore/scikit-bio/pull/1513))

### Bug fixes
* The `include_self` parameter was not being honored in `skbio.TreeNode.tips`. The scope of this bug was that if `TreeNode.tips` was called on a tip, it would always result in an empty `list` when unrolled.

* In `skbio.stats.ordination.ca`, `proportion_explained` was missing in the returned `OrdinationResults` object. ([#1345](https://github.com/biocore/scikit-bio/issues/1345))

* `skbio.diversity.beta_diversity` now handles qualitative metrics as expected such that `beta_diversity('jaccard', mat) == beta_diversity('jaccard', mat > 0)`. Please see [#1549](https://github.com/biocore/scikit-bio/issues/1549) for further detail.

* `skbio.stats.ordination.rda` The occasional column mismatch in output `biplot_scores` is fixed ([#1519](https://github.com/biocore/scikit-bio/issues/1519)).

### Deprecated functionality [stable]

### Deprecated functionality [experimental]

### Miscellaneous
* scikit-bio now depends on pandas >= 0.19.2, and is compatible with newer pandas versions (e.g. 0.20.3) that were previously incompatible.
* scikit-bio now depends on `numpy >= 1.9.2, < 1.14.0` for compatibility with Python 3.4, 3.5, and 3.6 and the available numpy conda packages in `defaults` and `conda-forge` channels.
* added support for running tests from `setup.py`. Both `python setup.py nosetests` and `python setup.py test` are now supported, however `python setup.py test` will only run a subset of the full test suite. ([#1341](https://github.com/biocore/scikit-bio/issues/1341))

## Version 0.5.1 (2016-11-12)

### Features
* Added `IntervalMetadata` and `Interval` classes in `skbio.metadata` to store, query, and manipulate information of a sub-region of a sequence. ([#1414](https://github.com/biocore/scikit-bio/issues/1414))
* `Sequence` and its child classes (including `GrammaredSequence`, `RNA`, `DNA`, `Protein`) now accept `IntervalMetadata` in their constructor API. Some of their relevant methods are also updated accordingly. ([#1430](https://github.com/biocore/scikit-bio/pull/1430))
* GenBank parser now reads and writes `Sequence` or its subclass objects with `IntervalMetadata`. ([#1440](https://github.com/biocore/scikit-bio/pull/1440))
* `DissimilarityMatrix` now has a new constructor method called `from_iterable`. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* `DissimilarityMatrix` now allows non-hollow matrices. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* `DistanceMatrix.from_iterable` now accepts a `validate=True` parameter. ([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* ``DistanceMatrix`` now has a new method called ``to_series`` to create a ``pandas.Series`` from a ``DistanceMatrix`` ([#1397](https://github.com/biocore/scikit-bio/issues/1397)).
* Added parallel beta diversity calculation support via `skbio.diversity.block_beta_diversity`. The issue and idea is discussed in ([#1181](https://github.com/biocore/scikit-bio/issues/1181), while the actual code changes are in [#1352](https://github.com/biocore/scikit-bio/pull/1352)).


### Backward-incompatible changes [stable]
* The constructor API for `Sequence` and its child classes (including `GrammaredSequence`, `RNA`, `DNA`, `Protein`) are changed from `(sequence, metadata=None, positional_metadata=None, lowercase=False)` to `(sequence, metadata=None, positional_metadata=None, interval_metadata=None, lowercase=False)`

  The changes are made to allow these classes to adopt `IntervalMetadata` object for interval features on the sequence. The `interval_metadata` parameter is added imediately after `positional_metadata` instead of appended to the end, because it is more natural and logical and, more importantly, because it is unlikely in practice to break user code. A user's code would break only if they had supplied `metadata`, `postional_metadata`, and `lowercase` parameters positionally. In the unlikely event that this happens, users will get an error telling them a bool isn't a valid `IntervalMetadata` type, so it won't silently produce buggy behavior.

### Backward-incompatible changes [experimental]
* Modifying basis handling in `skbio.stats.composition.ilr_inv` prior to checking for orthogonality.  Now the basis is strictly assumed to be in the Aitchison simplex.
* `DistanceMatrix.from_iterable` default behavior is now to validate matrix by computing all pairwise distances. Pass `validate=False` to get the previous behavior (no validation, but faster execution).([#1343](https://github.com/biocore/scikit-bio/issues/1343)).
* GenBank I/O now parses sequence features into the attribute of `interval_metadata` instead of `positiona_metadata`. And the key of `FEATURES` is removed from `metadata` attribute.

### Performance enhancements
* `TreeNode.shear` was rewritten for approximately a 25% performance increase. ([#1399](https://github.com/biocore/scikit-bio/pull/1399))
* The `IntervalMetadata` allows dramatic decrease in memory usage in reading GenBank files of feature rich sequences. ([#1159](https://github.com/biocore/scikit-bio/issues/1159))

### Bug fixes

* `skbio.tree.TreeNode.prune` and implicitly `skbio.tree.TreeNode.shear` were not handling a situation in which a parent was validly removed during pruning operations as may happen if the resulting subtree does not include the root. Previously, an `AttributeError` would raise as `parent` would be `None` in this situation.
* numpy linking was fixed for installation under El Capitan.
* A bug was introduced in #1398 into `TreeNode.prune` and fixed in #1416 in which, under the special case of a single descendent existing from the root, the resulting children parent references were not updated. The cause of the bug was a call made to `self.children.extend` as opposed to `self.extend` where the former is a `list.extend` without knowledge of the tree, while the latter is `TreeNode.extend` which is able to adjust references to `self.parent`.

### Miscellaneous

* Removed deprecated functions from `skbio.util`: `is_casava_v180_or_later`, `remove_files`, and `create_dir`.
* Removed deprecated `skbio.Sequence.copy` method.

## Version 0.5.0 (2016-06-14)

**IMPORTANT**: scikit-bio is no longer compatible with Python 2. scikit-bio is compatible with Python 3.4 and later.

### Features
* Added more descriptive error message to `skbio.io.registry` when attempting to read without specifying `into` and when there is no generator reader. ([#1326](https://github.com/biocore/scikit-bio/issues/1326))
* Added support for reference tags to `skbio.io.format.stockholm` reader and writer. ([#1348](https://github.com/biocore/scikit-bio/issues/1348))
* Expanded error message in `skbio.io.format.stockholm` reader when `constructor` is not passed, in order to provide better explanation to user. ([#1327](https://github.com/biocore/scikit-bio/issues/1327))
* Added `skbio.sequence.distance.kmer_distance` for computing the kmer distance between two sequences. ([#913](https://github.com/biocore/scikit-bio/issues/913))
* Added `skbio.sequence.Sequence.replace` for assigning a character to positions in a `Sequence`. ([#1222](https://github.com/biocore/scikit-bio/issues/1222))
* Added support for `pandas.RangeIndex`, lowering the memory footprint of default integer index objects. `Sequence.positional_metadata` and `TabularMSA.positional_metadata` now use `pd.RangeIndex` as the positional metadata index. `TabularMSA` now uses `pd.RangeIndex` as the default index. Usage of `pd.RangeIndex` over the previous `pd.Int64Index` [should be transparent](http://pandas.pydata.org/pandas-docs/version/0.18.0/whatsnew.html#range-index), so these changes should be non-breaking to users. scikit-bio now depends on pandas >= 0.18.0 ([#1308](https://github.com/biocore/scikit-bio/issues/1308))
* Added `reset_index=False` parameter to `TabularMSA.append` and `TabularMSA.extend` for resetting the MSA's index to the default index after appending/extending.
* Added support for partial pairwise calculations via `skbio.diversity.partial_beta_diversity`. ([#1221](https://github.com/biocore/scikit-bio/issues/1221), [#1337](https://github.com/biocore/scikit-bio/pull/1337)). This function is immediately deprecated as its return type will change in the future and should be used with caution in its present form (see the function's documentation for details).
* `TemporaryFile` and `NamedTemporaryFile` are now supported IO sources for `skbio.io` and related functionality.  ([#1291](https://github.com/biocore/scikit-bio/issues/1291))
* Added `tree_node_class=TreeNode` parameter to `skbio.tree.majority_rule` to support returning consensus trees of type `TreeNode` (the default) or a type that has the same interface as `TreeNode` (e.g. `TreeNode` subclasses) ([#1193](https://github.com/biocore/scikit-bio/pull/1193))
* `TreeNode.from_linkage_matrix` and `TreeNode.from_taxonomy` now support constructing `TreeNode` subclasses. `TreeNode.bifurcate` now supports `TreeNode` subclasses ([#1193](https://github.com/biocore/scikit-bio/pull/1193))
* The `ignore_metadata` keyword has been added to `TabularMSA.iter_positions` to improve performance when metadata is not necessary.
* Pairwise aligners in `skbio.alignment` now propagate per-sequence `metadata` objects (this does not include `positional_metadata`).

### Backward-incompatible changes [stable]

### Backward-incompatible changes [experimental]
* `TabularMSA.append` and `TabularMSA.extend` now require one of `minter`, `index`, or `reset_index` to be provided when incorporating new sequences into an MSA. Previous behavior was to auto-increment the index labels if `minter` and `index` weren't provided and the MSA had a default integer index, otherwise error. Use `reset_index=True` to obtain the previous behavior in a more explicit way.
* `skbio.stats.composition.ancom` now returns two `pd.DataFrame` objects, where it previously returned one. The first contains the ANCOM test results, as before, and the second contains percentile abundances of each feature in each group. The specific percentiles that are computed and returned is controlled by the new `percentiles` parameter to `skbio.stats.composition.ancom`. In the future, this second `pd.DataFrame` will not be returned by this function, but will be available through the [contingency table API](https://github.com/biocore/scikit-bio/issues/848). ([#1293](https://github.com/biocore/scikit-bio/issues/1293))
* `skbio.stats.composition.ancom` now performs multiple comparisons correction by default. The previous behavior of not performing multiple comparisons correction can be achieved by passing ``multiple_comparisons_correction=None``.
* The ``reject`` column in the first ``pd.DataFrame`` returned from `skbio.stats.composition.ancom` has been renamed ``Reject null hypothesis`` for clarity. ([#1375](https://github.com/biocore/scikit-bio/issues/1375))

### Bug fixes
* Fixed row and column names to `biplot_scores` in the `OrdinationResults` object from `skbio.stats.ordination`. This fix affect the `cca` and `rda` methods. ([#1322](https://github.com/biocore/scikit-bio/issues/1322))
* Fixed bug when using `skbio.io.format.stockholm` reader on file with multi-line tree with no id. Previously this raised an `AttributeError`, now it correctly handles this type of tree. ([#1334](https://github.com/biocore/scikit-bio/issues/1334))
* Fixed bug when reading Stockholm files with GF or GS features split over multiple lines. Previously, the feature text was simply concatenated because it was assumed to have trailing whitespace. There are examples of Stockholm files with and without trailing whitespace for multi-line features, so the `skbio.io.format.stockholm` reader now adds a single space when concatenating feature text without trailing whitespace to avoid joining words together. Multi-line trees stored as GF metadata are concatenated as they appear in the file; a space is not added when concatenating. ([#1328](https://github.com/biocore/scikit-bio/issues/1328))
* Fixed bug when using `Sequence.iter_kmers` on empty `Sequence` object. Previously this raised a `ValueError`, now it returns
an empty generator.
* Fixed minor bug where adding sequences to an empty `TabularMSA` with MSA-wide `positional_metadata` would result in a `TabularMSA` object in an inconsistent state. This could happen using `TabularMSA.append` or `TabularMSA.extend`. This bug only affects a `TabularMSA` object *without* sequences that has MSA-wide `positional_metadata` (for example, `TabularMSA([], positional_metadata={'column': []})`).
* `TreeNode.distance` now handles the situation in which `self` or `other` are ancestors. Previosly, a node further up the tree was used resulting in inflated distances. ([#807](https://github.com/biocore/scikit-bio/issues/807))
* `TreeNode.prune` can now handle a root with a single descendent. Previously, the root was ignored from possibly having a single descendent. ([#1247](https://github.com/biocore/scikit-bio/issues/1247))
* Providing the `format` keyword to `skbio.io.read` when creating a generator with an empty file will now return an empty generator instead of raising `StopIteration`. ([#1313](https://github.com/biocore/scikit-bio/issues/1313))
* `OrdinationResults` is now importable from `skbio` and `skbio.stats.ordination` and correctly linked from the documentation ([#1205](https://github.com/biocore/scikit-bio/issues/1205))
* Fixed performance bug in pairwise aligners resulting in 100x worse performance than in 0.2.4.

### Deprecated functionality [stable]
* Deprecated use of the term "non-degenerate", in favor of "definite". `GrammaredSequence.nondegenerate_chars`, `GrammaredSequence.nondegenerates`, and `GrammaredSequence.has_nondegenerates` have been renamed to `GrammaredSequence.definite_chars`, `GrammaredSequence.definites`, and `GrammaredSequence.has_definites`, respectively. The old names will be removed in scikit-bio 0.5.2. Relevant affected public classes include `GrammaredSequence`, `DNA`, `RNA`, and `Protein`.

### Deprecated functionality [experimental]
* Deprecated function `skbio.util.create_dir`. This function will be removed in scikit-bio 0.5.1. Please use the Python standard library
functionality described [here](https://docs.python.org/2/library/os.html#os.makedirs). ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated function `skbio.util.remove_files`. This function will be removed in scikit-bio 0.5.1. Please use the Python standard library
functionality described [here](https://docs.python.org/2/library/os.html#os.remove). ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated function `skbio.util.is_casava_v180_or_later`. This function will be removed in 0.5.1. Functionality moved to FASTQ sniffer.
([#833](https://github.com/biocore/scikit-bio/issues/833))

### Miscellaneous
* When installing scikit-bio via `pip`, numpy must now be installed first ([#1296](https://github.com/biocore/scikit-bio/issues/1296))

## Version 0.4.2 (2016-02-17)

Minor maintenance release. **This is the last Python 2.7 compatible release. Future scikit-bio releases will only support Python 3.**

### Features
* Added `skbio.tree.TreeNode.bifurcate` for converting multifurcating trees into bifurcating trees. ([#896](https://github.com/biocore/scikit-bio/issues/896))
* Added `skbio.io.format.stockholm` for reading Stockholm files into a `TabularMSA` and writing from a `TabularMSA`. ([#967](https://github.com/biocore/scikit-bio/issues/967))
* scikit-bio `Sequence` objects have better compatibility with numpy. For example, calling `np.asarray(sequence)` now converts the sequence to a numpy array of characters (the same as calling `sequence.values`).
* Added `skbio.sequence.distance` subpackage for computing distances between scikit-bio `Sequence` objects ([#913](https://github.com/biocore/scikit-bio/issues/913))
* Added ``skbio.sequence.GrammaredSequence``, which can be inherited from to create grammared sequences with custom alphabets (e.g., for use with TabularMSA) ([#1175](https://github.com/biocore/scikit-bio/issues/1175))
* Added ``skbio.util.classproperty`` decorator

### Backward-incompatible changes [stable]
* When sniffing or reading a file (`skbio.io.sniff`, `skbio.io.read`, or the object-oriented `.read()` interface), passing `newline` as a keyword argument to `skbio.io.open` now raises a `TypeError`. This backward-incompatible change to a stable API is necessary because it fixes a bug (more details in bug fix section below).
* When reading a FASTQ or QSEQ file and passing `variant='solexa'`, `ValueError` is now raised instead of `NotImplementedError`. This backward-incompatible change to a stable API is necessary to avoid creating a spin-locked process due to [a bug in Python](https://bugs.python.org/issue25786). See [#1256](https://github.com/biocore/scikit-bio/issues/1256) for details. This change is temporary and will be reverted to `NotImplementedError` when the bug is fixed in Python.

### Backward-incompatible changes [experimental]
* `skbio.io.format.genbank`: When reading GenBank files, the date field of the LOCUS line is no longer parsed into a `datetime.datetime` object and is left as a string. When writing GenBank files, the locus date metadata is expected to be a string instead of a `datetime.datetime` object ([#1153](https://github.com/biocore/scikit-bio/issues/1153))
* `Sequence.distance` now converts the input sequence (`other`) to its type before passing both sequences to `metric`. Previous behavior was to always convert to `Sequence`.

### Bug fixes
* Fixed bug when using `Sequence.distance` or `DistanceMatrix.from_iterable` to compute distances between `Sequence` objects with differing `metadata`/`positional_metadata` and passing `metric=scipy.spatial.distance.hamming` ([#1254](https://github.com/biocore/scikit-bio/issues/1254))
* Fixed performance bug when computing Hamming distances between `Sequence` objects in `DistanceMatrix.from_iterable` ([#1250](https://github.com/biocore/scikit-bio/issues/1250))
* Changed `skbio.stats.composition.multiplicative_replacement` to raise an error whenever a large value of `delta` is chosen ([#1241](https://github.com/biocore/scikit-bio/issues/1241))
* When sniffing or reading a file (`skbio.io.sniff`, `skbio.io.read`, or the object-oriented `.read()` interface), passing `newline` as a keyword argument to `skbio.io.open` now raises a `TypeError`. The file format's `newline` character will be used when opening the file. Previous behavior allowed overriding the format's `newline` character but this could cause issues with readers that assume newline characters are those defined by the file format (which is an entirely reasonable assumption). This bug is very unlikely to have surfaced in practice as the default `newline` behavior is *universal newlines mode*.
* DNA, RNA, and Protein are no longer inheritable because they assume an IUPAC alphabet.
* `DistanceMatrix` constructor provides more informative error message when data contains NaNs ([#1276](https://github.com/biocore/scikit-bio/issues/1276))

### Miscellaneous
* Warnings raised by scikit-bio now share a common subclass ``skbio.util.SkbioWarning``.

## Version 0.4.1 (2015-12-09)

### Features
* The ``TabularMSA`` object was added to represent and operate on tabular multiple sequence alignments. This satisfies [RFC 1](https://github.com/biocore/scikit-bio-rfcs/blob/master/active/001-tabular-msa.md). See the ``TabularMSA`` docs for full details.
* Added phylogenetic diversity metrics, including weighted UniFrac, unweighted UniFrac, and Faith's Phylogenetic Diversity. These are accessible as ``skbio.diversity.beta.unweighted_unifrac``, ``skbio.diversity.beta.weighted_unifrac``, and ``skbio.diversity.alpha.faith_pd``, respectively.
* Addition of the function ``skbio.diversity.alpha_diversity`` to support applying an alpha diversity metric to multiple samples in one call.
* Addition of the functions ``skbio.diversity.get_alpha_diversity_metrics`` and ``skbio.diversity.get_beta_diversity_metrics`` to support discovery of the alpha and beta diversity metrics implemented in scikit-bio.
* Added `skbio.stats.composition.ancom` function, a test for OTU differential abundance across sample categories. ([#1054](https://github.com/biocore/scikit-bio/issues/1054))
* Added `skbio.io.format.blast7` for reading BLAST+ output format 7 or BLAST output format 9 files into a `pd.DataFrame`. ([#1110](https://github.com/biocore/scikit-bio/issues/1110))
* Added `skbio.DissimilarityMatrix.to_data_frame` method for creating a ``pandas.DataFrame`` from a `DissimilarityMatrix` or `DistanceMatrix`. ([#757](https://github.com/biocore/scikit-bio/issues/757))
* Added support for one-dimensional vector of dissimilarities in `skbio.stats.distance.DissimilarityMatrix`
constructor. ([#6240](https://github.com/biocore/scikit-bio/issues/624))
* Added `skbio.io.format.blast6` for reading BLAST+ output format 6 or BLAST output format 8 files into a `pd.DataFrame`. ([#1110](https://github.com/biocore/scikit-bio/issues/1110))
* Added `inner`, `ilr`, `ilr_inv` and `clr_inv`, ``skbio.stats.composition``, which enables linear transformations on compositions ([#892](https://github.com/biocore/scikit-bio/issues/892)
* Added ``skbio.diversity.alpha.pielou_e`` function as an evenness metric of alpha diversity. ([#1068](https://github.com/biocore/scikit-bio/issues/1068))
* Added `to_regex` method to `skbio.sequence._iupac_sequence` ABC - it returns a regex object that matches all non-degenerate versions of the sequence.
* Added ``skbio.util.assert_ordination_results_equal`` function for comparing ``OrdinationResults`` objects in unit tests.
* Added ``skbio.io.format.genbank`` for reading and writing GenBank/GenPept for ``DNA``, ``RNA``, ``Protein`` and ``Sequence`` classes.
* Added ``skbio.util.RepresentationWarning`` for warning about substitutions, assumptions, or particular alterations that were made for the successful completion of a process.
* ``TreeNode.tip_tip_distances`` now supports nodes without an associated length. In this case, a length of 0.0 is assumed and an ``skbio.util.RepresentationWarning`` is raised. Previous behavior was to raise a ``NoLengthError``. ([#791](https://github.com/biocore/scikit-bio/issues/791))
* ``DistanceMatrix`` now has a new constructor method called `from_iterable`.
* ``Sequence`` now accepts ``lowercase`` keyword like ``DNA`` and others. Updated ``fasta``, ``fastq``, and ``qseq`` readers/writers for ``Sequence`` to reflect this.
* The ``lowercase`` method has been moved up to ``Sequence`` meaning all sequence objects now have a ``lowercase`` method.
* Added ``reverse_transcribe`` class method to ``RNA``.
* Added `Sequence.observed_chars` property for obtaining the set of observed characters in a sequence. ([#1075](https://github.com/biocore/scikit-bio/issues/1075))
* Added `Sequence.frequencies` method for computing character frequencies in a sequence. ([#1074](https://github.com/biocore/scikit-bio/issues/1074))
* Added experimental class-method ``Sequence.concat`` which will produce a new sequence from an iterable of existing sequences. Parameters control how positional metadata is propagated during a concatenation.
* ``TreeNode.to_array`` now supports replacing ``nan`` branch lengths in the resulting branch length vector with the value provided as ``nan_length_value``.
* ``skbio.io.format.phylip`` now supports sniffing and reading strict, sequential PHYLIP-formatted files into ``skbio.Alignment`` objects. ([#1006](https://github.com/biocore/scikit-bio/issues/1006))
* Added `default_gap_char` class property to ``DNA``, ``RNA``, and ``Protein`` for representing gap characters in a new sequence.

### Backward-incompatible changes [stable]
* `Sequence.kmer_frequencies` now returns a `dict`. Previous behavior was to return a `collections.Counter` if `relative=False` was passed, and a `collections.defaultdict` if `relative=True` was passed. In the case of a missing key, the `Counter` would return 0 and the `defaultdict` would return 0.0. Because the return type is now always a `dict`, attempting to access a missing key will raise a `KeyError`. This change *may* break backwards-compatibility depending on how the `Counter`/`defaultdict` is being used. We hope that in most cases this change will not break backwards-compatibility because both `Counter` and `defaultdict` are `dict` subclasses.

   If the previous behavior is desired, convert the `dict` into a `Counter`/`defaultdict`:

    ```python
    import collections
    from skbio import Sequence
    seq = Sequence('ACCGAGTTTAACCGAATA')

    # Counter
    freqs_dict = seq.kmer_frequencies(k=8)
    freqs_counter = collections.Counter(freqs_dict)

    # defaultdict
    freqs_dict = seq.kmer_frequencies(k=8, relative=True)
    freqs_default_dict = collections.defaultdict(float, freqs_dict)
    ```

   **Rationale:** We believe it is safer to return `dict` instead of `Counter`/`defaultdict` as this may prevent error-prone usage of the return value. Previous behavior allowed accessing missing kmers, returning 0 or 0.0 depending on the `relative` parameter. This is convenient in many cases but also potentially misleading. For example, consider the following code:

    ```python
    from skbio import Sequence
    seq = Sequence('ACCGAGTTTAACCGAATA')
    freqs = seq.kmer_frequencies(k=8)
    freqs['ACCGA']
    ```

    Previous behavior would return 0 because the kmer `'ACCGA'` is not present in the `Counter`. In one respect this is the correct answer because we asked for kmers of length 8; `'ACCGA'` is a different length so it is not included in the results. However, we believe it is safer to avoid this implicit behavior in case the user assumes there are no `'ACCGA'` kmers in the sequence (which there are!). A `KeyError` in this case is more explicit and forces the user to consider their query. Returning a `dict` will also be consistent with `Sequence.frequencies`.

### Backward-incompatible changes [experimental]
* Replaced ``PCoA``, ``CCA``, ``CA`` and ``RDA`` in ``skbio.stats.ordination`` with equivalent functions ``pcoa``, ``cca``, ``ca`` and ``rda``. These functions now take ``pd.DataFrame`` objects.
* Change ``OrdinationResults`` to have its attributes based on ``pd.DataFrame`` and ``pd.Series`` objects, instead of pairs of identifiers and values. The changes are as follows:
    - ``species`` and ``species_ids`` have been replaced by a ``pd.DataFrame`` named ``features``.
    - ``site`` and ``site_ids`` have been replaced by a ``pd.DataFrame`` named ``samples``.
    - ``eigvals`` is now a ``pd.Series`` object.
    - ``proportion_explained`` is now a ``pd.Series`` object.
    - ``biplot`` is now a ``pd.DataFrame`` object named ``biplot_scores``.
    - ``site_constraints`` is now a ``pd.DataFrame`` object named ``sample_constraints``.
* ``short_method_name`` and ``long_method_name`` are now required arguments of the ``OrdinationResults`` object.
* Removed `skbio.diversity.alpha.equitability`. Please use `skbio.diversity.alpha.pielou_e`, which is more accurately named and better documented. Note that `equitability` by default used logarithm base 2 while `pielou_e` uses logarithm base `e` as described in Heip 1974.
* ``skbio.diversity.beta.pw_distances`` is now called ``skbio.diversity.beta_diversity``. This function no longer defines a default metric, and ``metric`` is now the first argument to this function. This function can also now take a pairwise distances function as ``pairwise_func``.
* Deprecated function ``skbio.diversity.beta.pw_distances_from_table`` has been removed from scikit-bio as scheduled. Code that used this should be adapted to use ``skbio.diversity.beta_diversity``.
* ``TreeNode.index_tree`` now returns a 2-D numpy array as its second return value (the child node index) instead of a 1-D numpy array.
* Deprecated functions `skbio.draw.boxplots` and `skbio.draw.grouped_distributions` have been removed from scikit-bio as scheduled. These functions generated plots that were not specific to bioinformatics. These types of plots can be generated with seaborn or another general-purpose plotting package.
* Deprecated function `skbio.stats.power.bootstrap_power_curve` has been removed from scikit-bio as scheduled. Use `skbio.stats.power.subsample_power` or `skbio.stats.power.subsample_paired_power` followed by `skbio.stats.power.confidence_bound`.
* Deprecated function `skbio.stats.spatial.procrustes` has been removed from scikit-bio as scheduled in favor of `scipy.spatial.procrustes`.
* Deprecated class `skbio.tree.CompressedTrie` and function `skbio.tree.fasta_to_pairlist` have been removed from scikit-bio as scheduled in favor of existing general-purpose Python trie packages.
* Deprecated function `skbio.util.flatten` has been removed from scikit-bio as scheduled in favor of solutions available in the Python standard library (see [here](http://stackoverflow.com/a/952952/3639023) and [here](http://stackoverflow.com/a/406199/3639023) for examples).
* Pairwise alignment functions in `skbio.alignment` now return a tuple containing the `TabularMSA` alignment, alignment score, and start/end positions. The returned `TabularMSA`'s `index` is always the default integer index; sequence IDs are no longer propagated to the MSA. Additionally, the pairwise alignment functions now accept the following input types to align:
    - `local_pairwise_align_nucleotide`: `DNA` or `RNA`
    - `local_pairwise_align_protein`: `Protein`
    - `local_pairwise_align`: `IUPACSequence`
    - `global_pairwise_align_nucleotide`: `DNA`, `RNA`, or `TabularMSA[DNA|RNA]`
    - `global_pairwise_align_protein`: `Protein` or `TabularMSA[Protein]`
    - `global_pairwise_align`: `IUPACSequence` or `TabularMSA`
    - `local_pairwise_align_ssw`: `DNA`, `RNA`, or `Protein`. Additionally, this function now overrides the `protein` kwarg based on input type. `constructor` parameter was removed because the function now determines the return type based on input type.
* Removed `skbio.alignment.SequenceCollection` in favor of using a list or other standard library containers to store scikit-bio sequence objects (most `SequenceCollection` operations were simple list comprehensions). Use `DistanceMatrix.from_iterable` instead of `SequenceCollection.distances` (pass `key="id"` to exactly match original behavior).
* Removed `skbio.alignment.Alignment` in favor of `skbio.alignment.TabularMSA`.
* Removed `skbio.alignment.SequenceCollectionError` and `skbio.alignment.AlignmentError` exceptions as their corresponding classes no longer exist.

### Bug Fixes

* ``Sequence`` objects now handle slicing of empty positional metadata correctly. Any metadata that is empty will no longer be propagated by the internal ``_to`` constructor. ([#1133](https://github.com/biocore/scikit-bio/issues/1133))
* ``DissimilarityMatrix.plot()`` no longer leaves a white border around the
  heatmap it plots (PR #1070).
* TreeNode.root_at_midpoint`` no longer fails when a node with two equal length child branches exists in the tree. ([#1077](https://github.com/biocore/scikit-bio/issues/1077))
* ``TreeNode._set_max_distance``, as called through ``TreeNode.get_max_distance`` or ``TreeNode.root_at_midpoint`` would store distance information as ``list``s in the attribute ``MaxDistTips`` on each node in the tree, however, these distances were only valid for the node in which the call to ``_set_max_distance`` was made. The values contained in ``MaxDistTips`` are now correct across the tree following a call to ``get_max_distance``. The scope of impact of this bug is limited to users that were interacting directly with ``MaxDistTips`` on descendant nodes; this bug does not impact any known method within scikit-bio. ([#1223](https://github.com/biocore/scikit-bio/issues/1223))
* Added missing `nose` dependency to setup.py's `install_requires`. ([#1214](https://github.com/biocore/scikit-bio/issues/1214))
* Fixed issue that resulted in legends of ``OrdinationResult`` plots sometimes being truncated. ([#1210](https://github.com/biocore/scikit-bio/issues/1210))

### Deprecated functionality [stable]
* `skbio.Sequence.copy` has been deprecated in favor of `copy.copy(seq)` and `copy.deepcopy(seq)`.

### Miscellaneous
* Doctests are now written in Python 3.
* ``make test`` now validates MANIFEST.in using [check-manifest](https://github.com/mgedmin/check-manifest). ([#461](https://github.com/biocore/scikit-bio/issues/461))
* Many new alpha diversity equations added to ``skbio.diversity.alpha`` documentation. ([#321](https://github.com/biocore/scikit-bio/issues/321))
* Order of ``lowercase`` and ``validate`` keywords swapped in ``DNA``, ``RNA``, and ``Protein``.

## Version 0.4.0 (2015-07-08)

Initial beta release. In addition to the changes detailed below, the following
subpackages have been mostly or entirely rewritten and most of their APIs are
substantially different (and improved!):

* `skbio.sequence`
* `skbio.io`

The APIs of these subpackages are now stable, and all others are experimental. See the [API stability docs](https://github.com/biocore/scikit-bio/tree/0.4.0/doc/source/user/api_stability.rst) for more details, including what we mean by *stable* and *experimental* in this context. We recognize that this is a lot of backward-incompatible changes. To avoid these types of changes being a surprise to our users, our public APIs are now decorated to make it clear to developers when an API can be relied upon (stable) and when it may be subject to change (experimental).

### Features
* Added `skbio.stats.composition` for analyzing data made up of proportions
* Added new ``skbio.stats.evolve`` subpackage for evolutionary statistics. Currently contains a single function, ``hommola_cospeciation``, which implements a permutation-based test of correlation between two distance matrices.
* Added support for ``skbio.io.util.open_file`` and ``skbio.io.util.open_files`` to pull files from HTTP and HTTPS URLs. This behavior propagates to the I/O registry.
* FASTA/QUAL (``skbio.io.format.fasta``) and FASTQ (``skbio.io.format.fastq``) readers now allow blank or whitespace-only lines at the beginning of the file, between records, or at the end of the file. A blank or whitespace-only line in any other location will continue to raise an error [#781](https://github.com/biocore/scikit-bio/issues/781).
* scikit-bio now ignores leading and trailing whitespace characters on each line while reading FASTA/QUAL and FASTQ files.
* Added `ratio` parameter to `skbio.stats.power.subsample_power`. This allows the user to calculate power on groups for uneven size (For example, draw twice as many samples from Group B than Group A). If `ratio` is not set, group sizes will remain equal across all groups.
* Power calculations (`skbio.stats.power.subsample_power` and `skbio.stats.power.subsample_paired_power`) can use test functions that return multiple p values, like some multivariate linear regression models. Previously, the power calculations required the test to return a single p value.
* Added ``skbio.util.assert_data_frame_almost_equal`` function for comparing ``pd.DataFrame`` objects in unit tests.

### Performance enhancements
* The speed of quality score decoding has been significantly improved (~2x) when reading `fastq` files.
* The speed of `NucleotideSequence.reverse_complement` has been improved (~6x).

### Bug fixes
* Changed `Sequence.distance` to raise an error any time two sequences are passed of different lengths regardless of the `distance_fn` being passed. [(#514)](https://github.com/biocore/scikit-bio/issues/514)
* Fixed issue with ``TreeNode.extend`` where if given the children of another ``TreeNode`` object (``tree.children``), both trees would be left in an incorrect and unpredictable state. ([#889](https://github.com/biocore/scikit-bio/issues/889))
* Changed the way power was calculated in `subsample_paired_power` to move the subsample selection before the test is performed. This increases the number of Monte Carlo simulations performed during power estimation, and improves the accuracy of the returned estimate. Previous power estimates from `subsample_paired_power` should be disregarded and re-calculated. ([#910](https://github.com/biocore/scikit-bio/issues/910))
* Fixed issue where `randdm` was attempting to create asymmetric distance matrices.This was causing an error to be raised by the `DistanceMatrix` constructor inside of the `randdm` function, so that `randdm` would fail when attempting to create large distance matrices. ([#943](https://github.com/biocore/scikit-bio/issues/943))

### Deprecated functionality
* Deprecated `skbio.util.flatten`. This function will be removed in scikit-bio 0.3.1. Please use standard python library functionality
described here [Making a flat list out of lists of lists](http://stackoverflow.com/a/952952/3639023), [Flattening a shallow list](http://stackoverflow.com/a/406199/3639023) ([#833](https://github.com/biocore/scikit-bio/issues/833))
* Deprecated `skbio.stats.power.bootstrap_power_curve` will be removed in scikit-bio 0.4.1. It is deprecated in favor of using ``subsample_power`` or ``sample_paired_power`` to calculate a power matrix, and then the use of ``confidence_bounds`` to calculate the average and confidence intervals.

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
* `min_observations` was removed from `skbio.stats.power.subsample_power` and `skbio.stats.power.subsample_paired_power`. The minimum number of samples for subsampling depends on the data set and statistical tests. Having a default parameter to set unnecessary limitations on the technique.

### Miscellaneous
* Changed testing procedures
    - Developers should now use `make test`
    - Users can use `python -m skbio.test`
    - Added `skbio.util._testing.TestRunner` (available through `skbio.util.TestRunner`). Used to provide a `test` method for each module init file. This class represents a unified testing path which wraps all `skbio` testing functionality.
    - Autodetect Python version and disable doctests for Python 3.
* `numpy` is no longer required to be installed before installing scikit-bio!
* Upgraded checklist.py to check source files non-conforming to [new header style](http://scikit-bio.org/docs/latest/development/new_module.html). ([#855](https://github.com/biocore/scikit-bio/issues/855))
* Updated to use `natsort` >= 4.0.0.
* The method of subsampling was changed for ``skbio.stats.power.subsample_paired_power``. Rather than drawing a paired sample for the run and then subsampling for each count, the subsample is now drawn for each sample and each run. In test data, this did not significantly alter the power results.
* checklist.py now enforces `__future__` imports in .py files.

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
