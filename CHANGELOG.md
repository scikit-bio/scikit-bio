# scikit-bio changelog

## Version 0.1.4-dev (changes since 0.1.4 release go here)

### Features

* Added ability to compute distances between sequences in a ``SequenceCollection`` object ([#509](https://github.com/biocore/scikit-bio/issues/509)), and expanded ``Alignment.distance`` to allow the user to pass a function for computing distances (the default distance metric is still ``scipy.spatial.distance.hamming``) ([#194](https://github.com/biocore/scikit-bio/issues/194)).
* Added functionality to not penalize terminal gaps in global alignment. This functionality results in more biologically relevant global alignments (see [#537](https://github.com/biocore/scikit-bio/issues/537) for discussion of the issue) and is now the default behavior for global alignment.
* The python global aligners (``global_pairwise_align``, ``global_pairwise_align_nucleotide``, and ``global_pairwise_align_protein``) now support aligning pairs of sequences, pairs of alignments, and a sequence and an alignment (see [#550](https://github.com/biocore/scikit-bio/issues/550)). This functionality supports progressive multiple sequence alignment, among other things such as adding a sequence to an existing alignment. 

### Bug fixes

* Fixed bug that resulted in a ``ValueError`` from ``local_align_pairwise_nucleotide`` (see [#504](https://github.com/biocore/scikit-bio/issues/504)) under many circumstances. This would not generate incorrect results, but would cause the code to fail.

### Backward-incompatible changes

* Removed ``skbio.math``, leaving ``stats`` and ``diversity`` to become top level packages. For example, instead of ``from skbio.math.stats.ordination import PCoA`` you would now import ``from skbio.stats.ordination import PCoA``.
* The module ``skbio.math.gradient`` as well as the contents of ``skbio.math.subsample`` and ``skbio.math.stats.misc`` are now found in ``skbio.stats``. As an example, to import subsample: ``from skbio.stats import subsample``; to import everything from gradient: ``from skbio.stats.gradient import *``.
* The contents of ``skbio.math.stats.ordination.utils`` are now in ``skbio.stats.ordination``.
* Removed ``skbio.app`` subpackage (i.e., the *application controller framework*) as this code has been ported to the standalone [burrito](https://github.com/biocore/burrito) Python package. This code was not specific to bioinformatics and is useful for wrapping command-line applications in general.
* Removed ``skbio.core``, leaving ``alignment``, ``distance``, ``genetic_code``, ``sequence``, ``tree``, and ``workflow`` to become top level packages. For example, instead of ``from skbio.core.distance import DistanceMatrix`` you would now import ``from skbio.distance import DistanceMatrix``.
* Removed ``skbio.util.exception`` and ``skbio.util.warning`` (see [#577](https://github.com/biocore/scikit-bio/issues/577) for the reasoning behind this change). ``FileFormatError``, ``RecordError``, ``FieldError``, and ``EfficiencyWarning`` have been moved to ``skbio.util``. ``BiologicalSequenceError`` has been moved to ``skbio.sequence``. ``SequenceCollectionError`` and ``StockholmParseError`` have been moved to ``skbio.alignment``. ``DissimilarityMatrixError``, ``DistanceMatrixError``, ``DissimilarityMatrixFormatError``, and ``MissingIDError`` have been moved to ``skbio.distance``. ``TreeError``, ``NoLengthError``, ``DuplicateNodeError``, ``MissingNodeError``, and ``NoParentError`` have been moved to ``skbio.tree``. ``FastqParseError`` has been moved to ``skbio.parse.sequences``. ``GeneticCodeError``, ``GeneticCodeInitError``, and ``InvalidCodonError`` have been moved to ``skbio.genetic_code``.

### Miscellaneous

* Relaxed requirement in ``BiologicalSequence.distance`` that sequences being compared are of equal length. This is relevant for Hamming distance, so the check is still performed in that case, but other distance metrics may not have that requirement. See [#504](https://github.com/biocore/scikit-bio/issues/507)).
* Renamed ``powertrip.py`` repo-checking script to ``checklist.py`` for clarity.
* ``checklist.py`` now ensures that all unit tests import from a minimally deep API. For example, it will produce an error if ``skbio.core.distance.DistanceMatrix`` is used over ``skbio.DistanceMatrix``.

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
