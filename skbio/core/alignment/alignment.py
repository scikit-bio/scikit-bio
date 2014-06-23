#!/usr/bin/env python
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip, range
from future.utils import viewkeys

from collections import Counter, defaultdict
from warnings import warn

import numpy as np
from scipy.stats import entropy

from skbio.core.exception import SequenceCollectionError
from skbio.core.distance import DistanceMatrix


class SequenceCollection(object):
    """Class for storing collections of biological sequences.

    Parameters
    ----------
    seqs : list of `skbio.core.sequence.BiologicalSequence` objects
        The `skbio.core.sequence.BiologicalSequence` objects to load into
        a new `SequenceCollection` object.
    validate : bool, optional
        If True, runs the `is_valid` method after construction and raises
        `SequenceCollectionError` if ``is_valid == False``.

    Raises
    ------
    skbio.core.exception.SequenceCollectionError
        If ``validate == True`` and ``is_valid == False``.

    See Also
    --------
    skbio.core.sequence.BiologicalSequence
    skbio.core.sequence.NucelotideSequence
    skbio.core.sequence.DNASequence
    skbio.core.sequence.RNASequence
    Alignment
    skbio.parse.sequences
    skbio.parse.sequences.parse_fasta

    Examples
    --------
    >>> from skbio.core.alignment import SequenceCollection
    >>> from skbio.core.sequence import DNA
    >>> sequences = [DNA('ACCGT', id="seq1"),
    ...              DNA('AACCGGT', id="seq2")]
    >>> s1 = SequenceCollection(sequences)
    >>> s1
    <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

    """

    @classmethod
    def from_fasta_records(cls, fasta_records, seq_constructor,
                           validate=False):
        r"""Initialize a `SequenceCollection` object

        Parameters
        ----------
        fasta_records : iterator of tuples
            The records to load into a new `SequenceCollection` object. These
            should be tuples of ``(sequence_id, sequence)``.
        seq_constructor : skbio.core.sequence.BiologicalSequence
        validate : bool, optional
            If True, runs the `is_valid` method after construction and raises
            `SequenceCollectionError` if ``is_valid == False``.

        Returns
        -------
        SequenceCollection (or a derived class)
            The new `SequenceCollection` object.

        Raises
        ------
        skbio.core.exception.SequenceCollectionError
            If ``validate == True`` and ``is_valid == False``.

        See Also
        --------
        skbio.core.sequence.BiologicalSequence
        skbio.core.sequence.NucelotideSequence
        skbio.core.sequence.DNASequence
        skbio.core.sequence.RNASequence
        Alignment
        skbio.parse.sequences
        skbio.parse.sequences.parse_fasta

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.parse.sequences import parse_fasta
        >>> from StringIO import StringIO
        >>> from skbio.core.sequence import DNA
        >>> fasta_f = StringIO('>seq1\nACCGT\n>seq2\nAACCGGT\n')
        >>> s1 = SequenceCollection.from_fasta_records(
        ...     parse_fasta(fasta_f), DNA)
        >>> s1
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        >>> records = [('seq1', 'ACCGT'), ('seq2', 'AACCGGT')]
        >>> s1 = SequenceCollection.from_fasta_records(records, DNA)
        >>> s1
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        """
        data = []
        for seq_id, seq in fasta_records:
            try:
                id, description = seq_id.split(None, 1)
            except ValueError:
                id = seq_id.strip()
                description = None
            data.append(seq_constructor(seq, id=id,
                                        description=description))

        return cls(data, validate=validate)

    def __init__(self, seqs, validate=False):
        self._data = seqs
        self._id_to_index = {}
        for i, seq in enumerate(self._data):
            id = seq.id
            if id in self:
                raise SequenceCollectionError(
                    "All sequence ids must be unique, but "
                    "id %s is present multiple times." % id)
            else:
                self._id_to_index[seq.id] = i

        # This is bad because we're making a second pass through the sequence
        # collection to validate. We'll want to avoid this, but it's tricky
        # because different subclasses will want to define their own is_valid
        # methods.
        if validate and not self.is_valid():
            raise SequenceCollectionError(
                "%s failed to validate." % self.__class__.__name__)

    def __contains__(self, id):
        r"""The in operator.

        Parameters
        ----------
        id : str
            The id to look up in the `SequenceCollection`.

        Returns
        -------
        bool
            Indicates whether `id` corresponds to a sequence id
            in the `SequenceCollection`.

        .. shownumpydoc

        """
        return id in self._id_to_index

    def __eq__(self, other):
        r"""The equality operator.

        Parameters
        ----------
        other : `SequenceCollection`
            The `SequenceCollection` to test for equality against.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are equal.

        Notes
        -----
        `SequenceCollection` objects are equal if they are the same type,
        contain the same number of sequences, and if each of the
        `skbio.core.sequence.BiologicalSequence` objects, in order, are equal.

        .. shownumpydoc

        """
        if self.__class__ != other.__class__:
            return False
        elif len(self) != len(other):
            return False
        else:
            for self_seq, other_seq in zip(self, other):
                if self_seq != other_seq:
                    return False
        return True

    def __getitem__(self, index):
        r"""The indexing operator.

        Parameters
        ----------
        index : int, str
            The position or sequence id of the
            `skbio.core.sequence.BiologicalSequence` to return from the
            `SequenceCollection`.

        Returns
        -------
        `skbio.core.sequence.BiologicalSequence`
            The `skbio.core.sequence.BiologicalSequence` at the specified
            index in the `SequenceCollection`.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> s1[0]
        <DNASequence: ACCGT (length: 5)>
        >>> s1["seq1"]
        <DNASequence: ACCGT (length: 5)>

        .. shownumpydoc

        """
        if isinstance(index, str):
            return self.get_seq(index)
        else:
            return self._data[index]

    def __iter__(self):
        r"""The iter operator.

        Returns
        -------
        iterator
            `skbio.core.sequence.BiologicalSequence` iterator for the
            `SequenceCollection`.

        .. shownumpydoc

        """
        return iter(self._data)

    def __len__(self):
        r"""The len operator.

        Returns
        -------
        int
            The number of sequences in the `SequenceCollection`.

        .. shownumpydoc

        """
        return self.sequence_count()

    def __ne__(self, other):
        r"""The inequality operator.

        Parameters
        ----------
        other : `SequenceCollection`

        Returns
        -------
        bool
            Indicates whether self and other are not equal.

        Notes
        -----
        See `SequenceCollection.__eq__` for a description of what it means for
        a pair of `SequenceCollection` objects to be equal.

        .. shownumpydoc

        """
        return not self.__eq__(other)

    def __repr__(self):
        r"""The repr method.

        Returns
        -------
        str
            Returns a string representation of the object.

        Notes
        -----
        String representation contains the class name, the number of sequences
        in the `SequenceCollection` (n), and the mean and standard deviation
        sequence length.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(repr(s1))
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        .. shownumpydoc

        """
        cn = self.__class__.__name__
        count, center, spread = self.distribution_stats()
        return "<%s: n=%d; mean +/- std length=%.2f +/- %.2f>" \
            % (cn, count, center, spread)

    def __reversed__(self):
        """The reversed method.

        Returns
        -------
        iterator
            `skbio.core.sequence.BiologicalSequence` iterator for the
            `SequenceCollection` in reverse order.

        .. shownumpydoc

        """
        return reversed(self._data)

    def __str__(self):
        r"""The str method.

        Returns
        -------
        str
            Fasta-formatted string of all sequences in the object.

        .. shownumpydoc

        """
        return self.to_fasta()

    def distribution_stats(self, center_f=np.mean, spread_f=np.std):
        r"""Return sequence count, and center and spread of sequence lengths

        Parameters
        ----------
        center_f : function
            Should take a list-like object and return a single value
            representing the center of the distribution.
        spread_f : function
            Should take a list-like object and return a single value
            representing the spread of the distribution.

        Returns
        -------
        tuple of (int, float, float)
            The sequence count, center of length distribution, spread of length
            distribution.

        Notes
        -----
        Alternatives for `center_f` and `spread_f` could be median and median
        absolute deviation.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> s1.distribution_stats()
        (2, 6.0, 1.0)

        """
        if self.is_empty():
            return (0, 0.0, 0.0)
        else:
            sequence_count = self.sequence_count()
            sequence_lengths = self.sequence_lengths()
            return (sequence_count, center_f(sequence_lengths),
                    spread_f(sequence_lengths))

    def degap(self):
        r"""Return a new `SequenceCollection` with all gap characters removed.

        Returns
        -------
        SequenceCollection
            A new `SequenceCollection` where
            `skbio.core.sequence.BiologicalSequence.degap` has been called on
            each sequence.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> s2 = s1.degap()
        >>> s2
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

        """
        return SequenceCollection([seq.degap() for seq in self])

    def get_seq(self, id):
        r"""Return a sequence from the `SequenceCollection` by its id.

        Parameters
        ----------
        id, str
            The id of the sequence to return.

        Returns
        -------
        skbio.core.sequence.BiologicalSequence
            The `skbio.core.sequence.BiologicalSequence` with `id`.

        Raises
        ------
        KeyError
            If `id` is not in the `SequenceCollection` object.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1['seq1'])
        A--CCGT.

        """
        return self[self._id_to_index[id]]

    def ids(self):
        """Returns the `BiologicalSequence` ids

        Returns
        -------
        list
            The ordered list of ids for the
            `skbio.core.sequence.BiologicalSequence` objects in the
            `SequenceCollection`.

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.ids())
        ['seq1', 'seq2']

        """
        return [seq.id for seq in self]

    def int_map(self, prefix=""):
        """Create an integer-based mapping of sequence ids

        Parameters
        ----------
        prefix : str
            String prefix for new integer-based ids.

        Returns
        -------
        dict
            Mapping of new ids to sequences.
        dict
            Mapping of new ids to old ids.

        Notes
        -----
        This is useful when writing sequences out for use with programs that
        are picky about their sequence ids (e.g., raXML).

        The integer-based ids will be strings, for consistency (e.g., if prefix
        is passed) and begin at 1.

        References
        ----------
        RAxML Version 8: A tool for Phylogenetic Analysis and Post-Analysis of
        Large Phylogenies". In Bioinformatics, 2014

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> new_id_to_seqs, new_id_to_old_ids = s1.int_map()
        >>> print(repr(new_id_to_seqs['1']))
        <DNASequence: ACCGT (length: 5)>
        >>> print(repr(new_id_to_seqs['2']))
        <DNASequence: AACCGGT (length: 7)>
        >>> print(new_id_to_old_ids['1'])
        seq1
        >>> print(new_id_to_old_ids['2'])
        seq2

        """
        int_keys = []
        int_map = []
        for i, seq in enumerate(self):
            k = ("%s%d" % (prefix, i+1))
            int_map.append((k, seq))
            int_keys.append((k, seq.id))
        return dict(int_map), dict(int_keys)

    def is_empty(self):
        """Return True if the SequenceCollection is empty

        Returns
        -------
        bool
            ``True`` if `self` contains zero sequences, and ``False``
            otherwise.

        """
        return self.sequence_count() == 0

    def is_valid(self):
        """Return True if the SequenceCollection is valid

        Returns
        -------
        bool
            ``True`` if `self` is valid, and ``False`` otherwise.

        Notes
        -----
        Validity is defined as having no sequences containing characters
        outside of their valid character sets.

        See Also
        --------
        skbio.core.alignment.BiologicalSequence.is_valid

        Examples
        --------
        >>> from skbio.core.alignment import SequenceCollection
        >>> from skbio.core.sequence import DNA, RNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.is_valid())
        True
        >>> sequences = [RNA('ACCGT', id="seq1"),
        ...              RNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.is_valid())
        False

        """
        return self._validate_character_set()

    def iteritems(self):
        """Generator of id, sequence tuples

        Returns
        -------
        generator of tuples
            Each tuple contains ordered
            (`skbio.core.sequence.BiologicalSequence.id`,
            `skbio.core.sequence.BiologicalSequence`) pairs.

        """
        for seq in self:
            yield seq.id, seq

    def lower(self):
        """Converts all sequences to lowercase

        Returns
        -------
        SequenceCollection
            New `SequenceCollection` object where
            `skbio.core.sequence.BiologicalSequence.lower()` has been called
            on each sequence.

        See Also
        --------
        skbio.core.sequence.BiologicalSequence.lower
        upper

        """
        return self.__class__([seq.lower() for seq in self])

    def sequence_count(self):
        """Return the count of sequences in the `SequenceCollection`

        Returns
        -------
        int
            The number of sequences in the `SequenceCollection`.

        See Also
        --------
        sequence_lengths
        Alignment.sequence_length

        """
        return len(self._data)

    def k_word_frequencies(self, k, overlapping=True, constructor=str):
        """Return frequencies of length k words for sequences in Alignment

        Parameters
        ----------
        k : int
            The word length.
        overlapping : bool, optional
            Defines whether the k-words should be overlapping or not
            overlapping. This is only relevant when k > 1.
        constructor : type, optional
            The constructor for the returned k-words.

        Returns
        -------
        list
            List of ``collections.defaultdict`` objects, one for each sequence
            in the `Alignment`, representing the frequency of each character in
            each sequence of the `Alignment`.

        See Also
        --------
        position_frequencies

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('A', id="seq1"),
        ...              DNA('AT', id="seq2"),
        ...              DNA('TTTT', id="seq3")]
        >>> s1 = SequenceCollection(sequences)
        >>> for freqs in s1.k_word_frequencies(1):
        ...     print(freqs)
        defaultdict(<type 'int'>, {'A': 1.0})
        defaultdict(<type 'int'>, {'A': 0.5, 'T': 0.5})
        defaultdict(<type 'int'>, {'T': 1.0})
        >>> for freqs in s1.k_word_frequencies(2):
        ...     print(freqs)
        defaultdict(<type 'int'>, {})
        defaultdict(<type 'int'>, {'AT': 1.0})
        defaultdict(<type 'int'>, {'TT': 1.0})

        """
        result = []

        for s in self:
            result.append(s.k_word_frequencies(k, overlapping, constructor))
        return result

    def sequence_lengths(self):
        """Return lengths of the sequences in the `SequenceCollection`

        Returns
        -------
        list
            The ordered list of sequence lengths.

        See Also
        --------
        sequence_count

        """
        return [len(seq) for seq in self]

    def to_fasta(self):
        """Return fasta-formatted string representing the `SequenceCollection`

        Returns
        -------
        str
            A fasta-formatted string representing the `SequenceCollection`.

        See Also
        --------
        skbio.parse.sequences.parse_fasta
        """
        return ''.join([seq.to_fasta() for seq in self._data])

    def toFasta(self):
        """Return fasta-formatted string representing the `SequenceCollection`

        .. note:: Deprecated in skbio 0.3.0
                  `SequenceCollection.toFasta` will be removed in skbio 0.2.0,
                  it is replaced by `SequenceCollection.to_fasta` as the latter
                  adheres to PEP8 naming conventions. This is necessary to keep
                  in place now as these objects are sometimes passed into
                  code that expects a `cogent.alignment.Alignment` object
                  (e.g., PyNAST), so we need to support the method with this
                  name.

        Returns
        -------
        str
            A fasta-formatted string representing the `SequenceCollection`.

        """
        warn("SequenceCollection.toFasta() is deprecated. You should use "
             "SequenceCollection.to_fasta().")
        return self.to_fasta()

    def upper(self):
        """Converts all sequences to uppercase

        Returns
        -------
        SequenceCollection
            New `SequenceCollection` object where `BiologicalSequence.upper()`
            has been called on each sequence.

        See Also
        --------
        BiologicalSequence.upper
        lower

        """
        return self.__class__([seq.upper() for seq in self])

    def _validate_character_set(self):
        """Return ``True`` if all sequences are valid, ``False`` otherwise
        """
        for seq in self:
            if not seq.is_valid():
                return False
        return True


class Alignment(SequenceCollection):
    """Class for storing alignments of biological sequences.

    The ``Alignment`` class adds convenience methods to the
    ``SequenceCollection`` class to make it easy to work with alignments of
    biological sequences.

    Parameters
    ----------
    seqs : list of `skbio.core.sequence.BiologicalSequence` objects
        The `skbio.core.sequence.BiologicalSequence` objects to load into
        a new `Alignment` object.
    validate : bool, optional
        If True, runs the `is_valid` method after construction and raises
        `SequenceCollectionError` if ``is_valid == False``.
    score : float, optional
        The score of the alignment, if applicable (usually only if the
        alignment was just constructed).
    start_end_positions : iterable of two-item tuples, optional
        The start and end positions of each input sequence in the alignment,
        if applicable (usually only if the alignment was just constructed using
        a local alignment algorithm). Note that these should be indexes into
        the unaligned sequences, though the `Alignment` object itself doesn't
        know about these.

    Raises
    ------
    skbio.core.exception.SequenceCollectionError
        If ``validate == True`` and ``is_valid == False``.

    Notes
    -----
    By definition, all of the sequences in an alignment must be of the same
    length. For this reason, an alignment can be thought of as a matrix of
    sequences (rows) by positions (columns).

    See Also
    --------
    skbio.core.sequence.BiologicalSequence
    skbio.core.sequence.NucelotideSequence
    skbio.core.sequence.DNASequence
    skbio.core.sequence.RNASequence
    SequenceCollection
    skbio.parse.sequences
    skbio.parse.sequences.parse_fasta

    Examples
    --------
    >>> from skbio.core.alignment import Alignment
    >>> from skbio.core.sequence import DNA
    >>> sequences = [DNA('A--CCGT', id="seq1"),
    ...              DNA('AACCGGT', id="seq2")]
    >>> a1 = Alignment(sequences)
    >>> a1
    <Alignment: n=2; mean +/- std length=7.00 +/- 0.00>

    """

    def __init__(self, seqs, validate=False, score=None,
                 start_end_positions=None):
        super(Alignment, self).__init__(seqs, validate)

        if score is not None:
            self._score = float(score)
        self._start_end_positions = start_end_positions

    def distances(self):
        """Compute distances between all pairs of sequences

        Returns
        -------
        skbio.core.distance.DistanceMatrix
            Matrix containing the distances between all pairs of sequences.

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        skbio.core.distance.DistanceMatrix
        scipy.spatial.distance.hamming

        Notes
        -----
        Distances between sequences are computed as hamming distances, though
        this will be generalized (see #194).

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> seqs = [DNA("A-CCGGG", id="s1"),
        ...         DNA("ATCC--G", id="s2"),
        ...         DNA("ATCCGGA", id="s3")]
        >>> a1 = Alignment(seqs)
        >>> print(a1.distances())
        3x3 distance matrix
        IDs:
        s1, s2, s3
        Data:
        [[ 0.          0.42857143  0.28571429]
         [ 0.42857143  0.          0.42857143]
         [ 0.28571429  0.42857143  0.        ]]
        """
        sequence_count = self.sequence_count()
        dm = np.zeros((sequence_count, sequence_count))
        ids = []
        for i in range(sequence_count):
            self_i = self[i]
            ids.append(self_i.id)
            for j in range(i):
                dm[i, j] = dm[j, i] = self_i.distance(self[j])
        return DistanceMatrix(dm, ids)

    def score(self):
        """Returns the score of the alignment.

        Returns
        -------
        float, None
            The score of the alignment, or ``None`` if this was not provided on
            object construction.

        Notes
        -----
        This value will often be ``None``, as it is generally only going to be
        provided on construction if the alignment itself was built within
        scikit-bio.

        """
        return self._score

    def start_end_positions(self):
        """Returns the (start, end) positions for each aligned sequence.

        Returns
        -------
        list, None
            The ordered list of sequence lengths, or ``None`` if this was not
            provided on object construction.

        Notes
        -----
        The start/end positions indicate the range of the unaligned sequences
        that were aligned when the alignment was constructed. For example, if
        local alignment was performed on the sequences ACA and TACAT, depending
        on the specific algorithm that was used to perform the alignment the
        start/end positions would likely be: [(0,2), (1,3)], indicating that
        alignment looked like the following, and that the first and last
        positions of the second sequence were not included in the alignment:
        ACA
        ACA

        This value will often be ``None``, as it is generally only going to be
        provided on construction if the alignment itself was built within
        scikit-bio.

        """
        return self._start_end_positions

    def subalignment(self, seqs_to_keep=None, positions_to_keep=None,
                     invert_seqs_to_keep=False,
                     invert_positions_to_keep=False):
        """Returns new `Alignment` that is a subset of the current `Alignment`

        Parameters
        ----------
        seqs_to_keep : list, optional
            A list of sequence ids to be retained in the resulting
            `Alignment`. If this is not passed, the default will be to retain
            all sequences.
        positions_to_keep : list, optional
            A list of position ids to be retained in the resulting
            `Alignment`. If this is not passed, the default will be to retain
            all positions.
        invert_seqs_to_keep : bool, optional
            If `True`, the sequences identified in `seqs_to_keep` will be
            discarded, rather than retained.
        invert_positions_to_keep : bool, optional
            If `True`, the sequences identified in `positions_to_keep` will be
            discarded, rather than retained.

        Returns
        -------
        Alignment
            The specified subalignment.

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> seqs = [DNA("A-CCGGG", id="s1"),
        ...         DNA("ATCC--G", id="s2"),
        ...         DNA("ATCCGGA", id="s3")]
        >>> a1 = Alignment(seqs)
        >>> a1
        <Alignment: n=3; mean +/- std length=7.00 +/- 0.00>
        >>> a1.subalignment(seqs_to_keep=["s1", "s2"])
        <Alignment: n=2; mean +/- std length=7.00 +/- 0.00>
        >>> a1.subalignment(seqs_to_keep=["s1", "s2"],
        ...         invert_seqs_to_keep=True)
        <Alignment: n=1; mean +/- std length=7.00 +/- 0.00>
        >>> a1.subalignment(positions_to_keep=[0, 2, 3, 5])
        <Alignment: n=3; mean +/- std length=4.00 +/- 0.00>
        >>> a1.subalignment(positions_to_keep=[0, 2, 3, 5],
        ...         invert_positions_to_keep=True)
        <Alignment: n=3; mean +/- std length=3.00 +/- 0.00>
        >>> a1.subalignment(seqs_to_keep=["s1", "s2"],
        ...         positions_to_keep=[0, 2, 3, 5])
        <Alignment: n=2; mean +/- std length=4.00 +/- 0.00>

        """
        # if seqs_to_keep was not passed
        if seqs_to_keep is None:
            # and invert_seqs_to_keep is True
            if invert_seqs_to_keep:
                # return an empty alignment (because we're inverting the
                # default of keeping all sequences)
                return self.__class__([])
            # else if invert_seqs_to_keep is False
            else:
                # default to returning all sequences
                def keep_seq(i, id):
                    return True
        # else, if seqs_to_keep was passed
        else:
            seqs_to_keep = set(seqs_to_keep)
            # and invert_seqs_to_keep is True
            if invert_seqs_to_keep:
                # keep only sequences that were not listed in seqs_to_keep
                def keep_seq(i, id):
                    return not (id in seqs_to_keep or
                                i in seqs_to_keep)
            # else if invert_seqs_to_keep is False
            else:
                # keep only sequences that were listed in seqs_to_keep
                def keep_seq(i, id):
                    return (id in seqs_to_keep or
                            i in seqs_to_keep)

        # if positions_to_keep was not passed
        if positions_to_keep is None:
            # and invert_positions_to_keep is True
            if invert_positions_to_keep:
                # return an empty alignment (because we're inverting the
                # default of keeping all positions)
                return self.__class__([])
            # else if invert_positions_to_keep is False
            else:
                # default to returning all positions
                def keep_position(pos):
                    return True
        # else, if positions_to_keep was passed
        else:
            positions_to_keep = set(positions_to_keep)
            # and invert_positions_to_keep is True
            if invert_positions_to_keep:
                # keep only positions that were not listed in
                # positions_to_keep
                def keep_position(pos):
                    return pos not in positions_to_keep
            # else if invert_positions_to_keep is False
            else:
                # keep only sequences that were listed in positions_to_keep
                def keep_position(pos):
                    return pos in positions_to_keep

        # prep the result object
        result = []
        # iterate over sequences
        for sequence_index, seq in enumerate(self):
            # determine if we're keeping the current sequence
            if keep_seq(sequence_index, seq.id):
                # if so, iterate over the positions to determine which we're
                # keeping, and store them in a new list
                new_seq = [c for i, c in enumerate(seq) if keep_position(i)]
                # and then pack the resulting sequence into a new
                # BiologicalSequence object, of the same type as the current
                # object.
                # Note: This is bad, we are calling join too much. This
                # should be addressed in issue #194.
                result.append(seq.__class__(''.join(new_seq),
                              id=seq.id,
                              description=seq.description))
            # if we're not keeping the current sequence, move on to the next
            else:
                continue
        # pack the result up in the same type of object as the current object
        # and return it
        return self.__class__(result)

    def is_valid(self):
        """Return True if the Alignment is valid

        Returns
        -------
        bool
            ``True`` if `self` is valid, and ``False`` otherwise.

        Notes
        -----
        Validity is defined as having no sequences containing characters
        outside of their valid character sets, and all sequences being of equal
        length.

        See Also
        --------
        skbio.core.alignment.BiologicalSequence.is_valid

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA, RNA
        >>> sequences = [DNA('ACCGT--', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> a1 = Alignment(sequences)
        >>> a1.is_valid()
        True
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> a1 = Alignment(sequences)
        >>> print(a1.is_valid())
        False
        >>> sequences = [RNA('ACCGT--', id="seq1"),
        ...              RNA('AACCGGT', id="seq2")]
        >>> a1 = Alignment(sequences)
        >>> print(a1.is_valid())
        False

        """

        return super(Alignment, self).is_valid() and self._validate_lengths()

    def iter_positions(self, constructor=None):
        """Generator of Alignment positions (i.e., columns)

        Parameters
        ----------
        constructor : type, optional
            Constructor function for creating the positional values. By
            default, these will be the same type as corresponding
            `skbio.core.sequence.BiologicalSequence` in the
            `SequenceCollection` object, but you can pass a
            `skbio.core.sequence.BiologicalSequence` class here to ensure
            that they are all of consistent type, or ``str`` to have them
            returned as strings.

        Returns
        -------
        GeneratorType
            Generator of lists of positional values in the
            `SequenceCollection` (effectively the transpose of the alignment).

        See Also
        --------
        iter

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('ACCGT--', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> a1 = Alignment(sequences)
        >>> for position in a1.iter_positions():
        ...     print(position)
        [<DNASequence: A (length: 1)>, <DNASequence: A (length: 1)>]
        [<DNASequence: C (length: 1)>, <DNASequence: A (length: 1)>]
        [<DNASequence: C (length: 1)>, <DNASequence: C (length: 1)>]
        [<DNASequence: G (length: 1)>, <DNASequence: C (length: 1)>]
        [<DNASequence: T (length: 1)>, <DNASequence: G (length: 1)>]
        [<DNASequence: - (length: 1)>, <DNASequence: G (length: 1)>]
        [<DNASequence: - (length: 1)>, <DNASequence: T (length: 1)>]

        >>> for position in a1.iter_positions(constructor=str):
        ...     print(position)
        ['A', 'A']
        ['C', 'A']
        ['C', 'C']
        ['G', 'C']
        ['T', 'G']
        ['-', 'G']
        ['-', 'T']

        """
        if constructor is None:
            def constructor(s):
                return s
        for i in range(self.sequence_length()):
            position = [constructor(seq[i]) for seq in self]
            yield position

    def majority_consensus(self, constructor=None):
        """Return the majority consensus sequence for the `Alignment`

        Parameters
        ----------
        constructor : function, optional
            Constructor function for creating the consensus sequence. By
            default, this will be the same type as the first sequence in the
            `Alignment`.

        Returns
        -------
        skbio.core.sequence.BiologicalSequence
            The consensus sequence of the `Alignment`. In other words, at each
            position the most common character is chosen, and those characters
            are combined to create a new sequence.

        Notes
        -----
        If there are two characters that are equally abundant in the sequence
        at a given position, the choice of which of those characters will be
        present at that position in the result is arbitrary.

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> a1.majority_consensus()
        <DNASequence: AT-C (length: 4)>
        >>> a1.majority_consensus(constructor=str)
        'AT-C'

        """
        # handle empty Alignment case
        if self.is_empty():
            return ''

        if constructor is None:
            constructor = self[0].__class__
        result = []
        for c in self.position_counters():
            # Counter.most_common returns an ordered list of the
            # n most common (sequence, count) items in Counter. Here
            # we set n=1, and take only the character, not the count.
            result.append(c.most_common(1)[0][0])
        result = ''.join(result)
        return constructor(result)

    def omit_gap_positions(self, maximum_gap_frequency):
        """Returns Alignment with positions filtered based on gap frequency

        Parameters
        ----------
        maximum_gap_frequency : float
            The maximum fraction of the sequences that can contain a gap at a
            given position for that position to be retained in the resulting
            `Alignment`.

        Returns
        -------
        Alignment
            The subalignment containing only the positions with gaps in fewer
            than `maximum_gap_frequency` fraction of the sequences.

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> a2 = a1.omit_gap_positions(0.50)
        >>> a2
        <Alignment: n=3; mean +/- std length=3.00 +/- 0.00>
        >>> print(a2[0])
        AC-
        >>> print(a2[1])
        ATC
        >>> print(a2[2])
        TTC

        """
        # handle empty Alignment case
        if self.is_empty():
            return self.__class__([])

        position_frequencies = self.position_frequencies()
        gap_alphabet = self[0].gap_alphabet()

        positions_to_keep = []
        for i, f in enumerate(position_frequencies):
            gap_frequency = sum([f[c] for c in gap_alphabet])
            if gap_frequency <= maximum_gap_frequency:
                positions_to_keep.append(i)
        return self.subalignment(positions_to_keep=positions_to_keep)

    def omit_gap_sequences(self, maximum_gap_frequency):
        """Returns Alignment with sequences filtered based on gap frequency

        Parameters
        ----------
        maximum_gap_frequency : float
            The maximum fraction of the positions that can contain a gap in a
            given sequence for that sequence to be retained in the resulting
            `Alignment`.

        Returns
        -------
        Alignment
            The subalignment containing only the sequences with gaps in fewer
            than `maximum_gap_frequency` fraction of the positions.

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> a2 = a1.omit_gap_sequences(0.49)
        >>> a2
        <Alignment: n=2; mean +/- std length=4.00 +/- 0.00>
        >>> print(a2[0])
        AT-C
        >>> print(a2[1])
        TT-C

        """
        # handle empty Alignment case
        if self.is_empty():
            return self.__class__([])

        base_frequencies = self.k_word_frequencies(k=1)
        gap_alphabet = self[0].gap_alphabet()
        seqs_to_keep = []
        for seq, f in zip(self, base_frequencies):
            gap_frequency = sum([f[c] for c in gap_alphabet])
            if gap_frequency <= maximum_gap_frequency:
                seqs_to_keep.append(seq.id)
        return self.subalignment(seqs_to_keep=seqs_to_keep)

    def position_counters(self):
        """Return collection.Counter object for positions in Alignment

        Returns
        -------
        list
            List of ``collection.Counter`` objects, one for each position in
            the `Alignment`.

        See Also
        --------
        position_frequencies
        position_entropies

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> for counter in a1.position_counters():
        ...     print(counter)
        Counter({'A': 2, 'T': 1})
        Counter({'T': 2, 'C': 1})
        Counter({'-': 3})
        Counter({'C': 2, '-': 1})

        """
        return [Counter(p) for p in self.iter_positions(constructor=str)]

    def position_frequencies(self):
        """Return frequencies of characters for positions in Alignment

        Returns
        -------
        list
            List of ``collection.defaultdict`` objects, one for each position
            in the `Alignment`, representing the frequency of each character in
            the `Alignment` at that position.

        See Also
        --------
        position_counters
        position_entropies
        k_word_frequencies

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> position_freqs = a1.position_frequencies()
        >>> print(round(position_freqs[0]['A'],3))
        0.667
        >>> print(round(position_freqs[1]['A'],3))
        0.0

        """
        result = []
        # handle the empty Alignment case
        if self.is_empty():
            return result

        count = 1 / self.sequence_count()
        for p in self.iter_positions(constructor=str):
            current_freqs = defaultdict(float)
            for c in p:
                current_freqs[c] += count
            result.append(current_freqs)
        return result

    def position_entropies(self, base=None,
                           nan_on_non_standard_chars=True):
        """Return Shannon entropy of positions in Alignment

        Parameters
        ----------
        base : float, optional
            log base for entropy calculation. If not passed, default will be e
            (i.e., natural log will be computed).
        nan_on_non_standard_chars : bool, optional
            if True, the entropy at positions containing characters outside of
            the first sequence's `iupac_standard_characters` will be `np.nan`.
            This is useful, and the default behavior, as it's not clear how a
            gap or degenerate character should contribute to a positional
            entropy. This issue was described in [1]_.

        Returns
        -------
        list
            List of floats of Shannon entropy at `Alignment` positions. Shannon
            entropy is defined in [2]_.

        See Also
        --------
        position_counters
        position_frequencies

        References
        ----------
        .. [1] Identifying DNA and protein patterns with statistically
           significant alignments of multiple sequences.
           Hertz GZ, Stormo GD.
           Bioinformatics. 1999 Jul-Aug;15(7-8):563-77.
        .. [2] A Mathematical Theory of Communication
           CE Shannon
           The Bell System Technical Journal (1948).

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> print(a1.position_entropies())
        [0.63651416829481278, 0.63651416829481278, nan, nan]

        """
        result = []
        # handle empty Alignment case
        if self.is_empty():
            return result

        iupac_standard_characters = self[0].iupac_standard_characters()
        for f in self.position_frequencies():
            if (nan_on_non_standard_chars and
                    len(viewkeys(f) - iupac_standard_characters) > 0):
                result.append(np.nan)
            else:
                result.append(entropy(list(f.values()), base=base))
        return result

    def sequence_length(self):
        """Return the number of positions in Alignment

        Returns
        -------
        int
            The number of positions in `Alignment`.

        See Also
        --------
        sequence_lengths
        sequence_count

        Examples
        --------
        >>> from skbio.core.alignment import Alignment
        >>> from skbio.core.sequence import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> a1.sequence_length()
        4

        """
        # handle the empty Alignment case
        if self.is_empty():
            return 0
        else:
            return len(self._data[0])

    def to_phylip(self, map_labels=False, label_prefix=""):
        """Return phylip-formatted string representing the `SequenceCollection`

        Returns
        -------
        str
            A phylip-formatted string representing the `SequenceCollection`.

        """
        if not self._validate_lengths():
            raise SequenceCollectionError("PHYLIP-formatted string can only "
                                          "be generated if all sequences are "
                                          "of equal length.")

        if self.is_empty():
            raise SequenceCollectionError("PHYLIP-formatted string can only "
                                          "be generated if there is at least "
                                          "one sequence in the Alignment.")

        sequence_length = self.sequence_length()
        if sequence_length == 0:
            raise SequenceCollectionError("PHYLIP-formatted string can only "
                                          "be generated if there is at least "
                                          "one position in the Alignment.")

        ids = self.ids()
        sequence_count = self.sequence_count()
        result = ["%d %d" % (sequence_count, sequence_length)]
        if map_labels:
            _, new_id_to_old_id = self.int_map(prefix=label_prefix)
            old_id_to_new_id = {v: k for k, v in new_id_to_old_id.items()}
        else:
            new_id_to_old_id = {seq_id: seq_id for seq_id in ids}
            old_id_to_new_id = new_id_to_old_id

        for seq_id in ids:
            new_id = old_id_to_new_id[seq_id]
            seq = self[seq_id]
            result.append("%s %s" % (new_id, str(seq)))

        return '\n'.join(result), new_id_to_old_id

    def _validate_lengths(self):
        """Return ``True`` if all sequences same length, ``False`` otherwise
        """
        seq1_length = self.sequence_length()
        for seq in self:
            if seq1_length != len(seq):
                return False
        return True
