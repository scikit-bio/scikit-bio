# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip, range
from future.utils import viewkeys, viewitems
from six import StringIO

import operator
from collections import Counter, defaultdict, OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import entropy

from skbio.alignment._text import add_letter
from skbio._base import SkbioObject
from skbio.sequence import Sequence
from skbio.stats.distance import DistanceMatrix
from skbio.io.util import open_file
from ._exception import (SequenceCollectionError, StockholmParseError,
                         AlignmentError)


class SequenceCollection(SkbioObject):
    """Class for storing collections of biological sequences.

    Parameters
    ----------
    seqs : list of `skbio.Sequence` objects
        The `skbio.Sequence` objects to load into a new `SequenceCollection`
        object.
    validate : bool, optional
        If True, runs the `is_valid` method after construction and raises
        `SequenceCollectionError` if ``is_valid == False``.

    Raises
    ------
    skbio.SequenceCollectionError
        If ``validate == True`` and ``is_valid == False``.

    See Also
    --------
    skbio
    skbio.DNA
    skbio.RNA
    skbio.Protein
    Alignment

    Examples
    --------
    >>> from skbio import SequenceCollection
    >>> from skbio import DNA
    >>> sequences = [DNA('ACCGT', id="seq1"),
    ...              DNA('AACCGGT', id="seq2")]
    >>> s1 = SequenceCollection(sequences)
    >>> s1
    <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

    """
    default_write_format = 'fasta'

    def __init__(self, seqs):
        # TODO: find a good way to support generic Sequence objects in
        # SequenceCollection and Alignment. The issue is that some methods
        # assume that a sequence has knowledge of gap characters and a
        # standard alphabet, which aren't present on Sequence. For now, if
        # these methods are called by a user they'll get an error (likely
        # an AttributeError).
        self._data = seqs
        self._id_to_index = {}
        for i, seq in enumerate(self._data):
            id = seq.id
            if id in self:
                raise SequenceCollectionError(
                    "All sequence ids must be unique, but "
                    "id '%s' is present multiple times." % id)
            else:
                self._id_to_index[seq.id] = i

    def __contains__(self, id):
        r"""The in operator.

        Parameters
        ----------
        id : str
            The `skbio.Sequence.id` to look up in the `SequenceCollection`.

        Returns
        -------
        bool
            Returns `True` if `id` is the `skbio.Sequence.id` of a sequence in
            the `SequenceCollection`.

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
        `skbio.Sequence` objects, in order, are equal.

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
            The position or sequence id of the `skbio.Sequence` to return from
            the `SequenceCollection`.

        Returns
        -------
        skbio.Sequence
            The `skbio.Sequence` at the specified index in the
            `SequenceCollection`.

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> s1[0]
        DNA('ACCGT', length=5, id='seq1')
        >>> s1["seq1"]
        DNA('ACCGT', length=5, id='seq1')

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
            `skbio.Sequence` iterator for the `SequenceCollection`.

        """
        return iter(self._data)

    def __len__(self):
        r"""The len operator.

        Returns
        -------
        int
            The number of sequences in the `SequenceCollection`.

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
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(repr(s1))
        <SequenceCollection: n=2; mean +/- std length=6.00 +/- 1.00>

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
            `skbio.Sequence` iterator for the `SequenceCollection` in reverse
            order.

        """
        return reversed(self._data)

    def __str__(self):
        r"""The str method.

        Returns
        -------
        str
            Fasta-formatted string of all sequences in the object.

        """
        fh = StringIO()
        self.write(fh, format='fasta')
        fasta_str = fh.getvalue()
        fh.close()
        return fasta_str

    def distances(self, distance_fn):
        """Compute distances between all pairs of sequences

        Parameters
        ----------
        distance_fn : function
            Function for computing the distance between a pair of sequences.
            This must take two sequences as input (as `skbio.Sequence` objects)
            and return a single integer or float value.

        Returns
        -------
        skbio.DistanceMatrix
            Matrix containing the distances between all pairs of sequences.

        See Also
        --------
        skbio.DistanceMatrix
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from scipy.spatial.distance import hamming
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> seqs = [DNA("ACCGGGTT", id="s1"),
        ...         DNA("ACTTGGTT", id="s2"),
        ...         DNA("ACTAGGTT", id="s3")]
        >>> a1 = SequenceCollection(seqs)
        >>> print(a1.distances(hamming))
        3x3 distance matrix
        IDs:
        's1', 's2', 's3'
        Data:
        [[ 0.     0.25   0.25 ]
         [ 0.25   0.     0.125]
         [ 0.25   0.125  0.   ]]

        """
        sequence_count = self.sequence_count()
        dm = np.zeros((sequence_count, sequence_count))
        ids = []
        for i in range(sequence_count):
            self_i = self[i]
            ids.append(self_i.id)
            for j in range(i):
                dm[i, j] = dm[j, i] = self_i.distance(self[j], distance_fn)
        return DistanceMatrix(dm, ids)

    def distribution_stats(self, center_f=np.mean, spread_f=np.std):
        r"""Return sequence count, and center and spread of sequence lengths

        Parameters
        ----------
        center_f : function
            Should take an array_like object and return a single value
            representing the center of the distribution.
        spread_f : function
            Should take an array_like object and return a single value
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
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
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
            A new `SequenceCollection` where `skbio.Sequence.degap` has been
            called on each sequence.

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
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
        id : str
            The id of the sequence to return.

        Returns
        -------
        skbio.Sequence
            The `skbio.Sequence` with `id`.

        Raises
        ------
        KeyError
            If `id` is not in the `SequenceCollection` object.

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1['seq1'])
        A--CCGT.

        """
        return self[self._id_to_index[id]]

    def ids(self):
        """Returns the `Sequence` ids

        Returns
        -------
        list
            The ordered list of ids for the `skbio.Sequence` objects in the
            `SequenceCollection`.

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.ids())
        ['seq1', 'seq2']

        """
        return [seq.id for seq in self]

    def update_ids(self, ids=None, func=None, prefix=""):
        """Update sequence IDs on the sequence collection.

        IDs can be updated by providing a sequence of new IDs (`ids`) or a
        function that maps current IDs to new IDs (`func`).

        Default behavior (if `ids` and `func` are not provided) is to create
        new IDs that are unique postive integers (starting at 1) cast as
        strings, optionally preceded by `prefix`. For example, ``('1', '2',
        '3', ...)``.

        Parameters
        ----------
        ids : sequence of str, optional
            New IDs to update on the sequence collection.
        func : function, optional
            Function accepting a sequence of current IDs and returning a
            sequence of new IDs to update on the sequence collection.
        prefix : str, optional
            If `ids` and `func` are both ``None``, `prefix` is prepended to
            each new integer-based ID (see description of default behavior
            above).

        Returns
        -------
        SequenceCollection
            New ``SequenceCollection`` (or subclass) containing sequences with
            updated IDs.
        dict
            Mapping of new IDs to old IDs.

        Raises
        ------
        SequenceCollectionError
            If both `ids` and `func` are provided, `prefix` is provided with
            either `ids` or `func`, or the number of new IDs does not match the
            number of sequences in the sequence collection.

        Notes
        -----
        The default behavior can be useful when writing sequences out for use
        with programs that are picky about their sequence IDs
        (e.g., RAxML [1]_).

        References
        ----------
        .. [1] RAxML Version 8: A tool for Phylogenetic Analysis and
           Post-Analysis of Large Phylogenies". In Bioinformatics, 2014

        Examples
        --------
        Define a sequence collection containing two sequences with IDs "abc"
        and "def":

        >>> from skbio import DNA, SequenceCollection
        >>> sequences = [DNA('A--CCGT.', id="abc"),
        ...              DNA('.AACCG-GT.', id="def")]
        >>> s1 = SequenceCollection(sequences)
        >>> s1.ids()
        ['abc', 'def']

        Update the IDs in the sequence collection, obtaining a new sequence
        collection with IDs that are integer-based:

        >>> s2, new_to_old_ids = s1.update_ids()
        >>> s2.ids()
        ['1', '2']

        Alternatively, we can specify a function to map the current IDs to new
        IDs. Let's define a function that appends ``'-new'`` to each ID:

        >>> def id_mapper(ids):
        ...     return [id_ + '-new' for id_ in ids]
        >>> s3, new_to_old_ids = s1.update_ids(func=id_mapper)
        >>> s3.ids()
        ['abc-new', 'def-new']

        We can also directly update the IDs with a new sequence of IDs:

        >>> s4, new_to_old_ids = s1.update_ids(ids=['ghi', 'jkl'])
        >>> s4.ids()
        ['ghi', 'jkl']

        """
        if ids is not None and func is not None:
            raise SequenceCollectionError("ids and func cannot both be "
                                          "provided.")
        if (ids is not None and prefix) or (func is not None and prefix):
            raise SequenceCollectionError("prefix cannot be provided if ids "
                                          "or func is provided.")

        if ids is not None:
            def func(_):
                return ids

        elif func is None:
            def func(_):
                new_ids = []
                for i in range(1, len(self) + 1):
                    new_ids.append("%s%d" % (prefix, i))
                return new_ids

        old_ids = self.ids()
        new_ids = func(old_ids)

        if len(new_ids) != len(old_ids):
            raise SequenceCollectionError(
                "Number of new IDs must be equal to the number of existing "
                "IDs (%d != %d)." % (len(new_ids), len(old_ids)))

        new_to_old_ids = dict(zip(new_ids, old_ids))

        new_seqs = []
        for new_id, seq in zip(new_ids, self):
            # HACK: Sequence objects are currently immutable. We used to have a
            # Sequence.copy/to method that created a new Sequence object,
            # optionally with attribute(s) set to new values. Here we want to
            # retain all state in the Sequence object and only update the ID.
            # There is a private method `_to` that accomplishes this for us.
            # This method used to be public and is fully tested and documented.
            # In the future, Sequence objects will have a `copy` method and
            # their attributes will be reassignable. When that happens, this
            # code can be updated to something like:
            #
            #     new_seq = seq.copy()
            #     new_seq.id = new_id
            #     new_seqs.append(new_seq)
            new_seqs.append(seq._to(id=new_id))

        return self.__class__(new_seqs), new_to_old_ids

    def is_empty(self):
        """Return True if the SequenceCollection is empty

        Returns
        -------
        bool
            ``True`` if `self` contains zero sequences, and ``False``
            otherwise.

        """
        return self.sequence_count() == 0

    def iteritems(self):
        """Generator of id, sequence tuples

        Returns
        -------
        generator of tuples
            Each tuple contains ordered (`skbio.Sequence.id`, `skbio.Sequence`)
            pairs.

        """
        for seq in self:
            yield seq.id, seq

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

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('A--CCGT.', id="seq1"),
        ...              DNA('.AACCG-GT.', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.sequence_count())
        2

        """
        return len(self._data)

    def kmer_frequencies(self, k, overlap=True, relative=False):
        """Return k-word frequencies for sequences in ``SequenceCollection``.

        Parameters
        ----------
        k : int
            The word length.
        overlapping : bool, optional
            Defines whether the k-words should be overlapping or not
            overlapping. This is only relevant when `k` > 1.

        Returns
        -------
        list
            List of ``collections.Counter`` objects, one for each sequence
            in the ``SequenceCollection``, representing the frequency of each
            k-word in each sequence of the ``SequenceCollection``.

        See Also
        --------
        Alignment.position_frequencies

        Examples
        --------
        >>> from skbio import SequenceCollection, DNA
        >>> sequences = [DNA('A', id="seq1"),
        ...              DNA('AT', id="seq2"),
        ...              DNA('TTTT', id="seq3")]
        >>> s1 = SequenceCollection(sequences)
        >>> for freqs in s1.kmer_frequencies(1):
        ...     print(freqs)
        Counter({'A': 1})
        Counter({'A': 1, 'T': 1})
        Counter({'T': 4})
        >>> for freqs in s1.kmer_frequencies(2):
        ...     print(freqs)
        Counter()
        Counter({'AT': 1})
        Counter({'TT': 3})

        """
        return [s.kmer_frequencies(k, overlap=overlap, relative=relative)
                for s in self]

    def sequence_lengths(self):
        """Return lengths of the sequences in the `SequenceCollection`

        Returns
        -------
        list
            The ordered list of sequence lengths.

        See Also
        --------
        sequence_count

        Examples
        --------
        >>> from skbio import SequenceCollection
        >>> from skbio import DNA
        >>> sequences = [DNA('ACCGT', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> s1 = SequenceCollection(sequences)
        >>> print(s1.sequence_lengths())
        [5, 7]

        """
        return [len(seq) for seq in self]


class Alignment(SequenceCollection):
    """Class for storing alignments of biological sequences.

    The ``Alignment`` class adds methods to the ``SequenceCollection`` class
    that are useful for working with aligned biological sequences.

    Parameters
    ----------
    seqs : list of `skbio.Sequence` objects
        The `skbio.Sequence` objects to load into a new `Alignment` object.
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
        know about these unless it is degapped.

    Raises
    ------
    skbio.SequenceCollectionError
        If ``validate == True`` and ``is_valid == False``.
    skbio.AlignmentError
        If not all the sequences have the same length.

    Notes
    -----
    By definition, all of the sequences in an alignment must be of the same
    length. For this reason, an alignment can be thought of as a matrix of
    sequences (rows) by positions (columns).

    See Also
    --------
    skbio
    skbio.DNA
    skbio.RNA
    skbio.Protein
    SequenceCollection

    Examples
    --------
    >>> from skbio import Alignment
    >>> from skbio import DNA
    >>> sequences = [DNA('A--CCGT', id="seq1"),
    ...              DNA('AACCGGT', id="seq2")]
    >>> a1 = Alignment(sequences)
    >>> a1
    <Alignment: n=2; mean +/- std length=7.00 +/- 0.00>

    """

    def __init__(self, seqs, score=None, start_end_positions=None):
        super(Alignment, self).__init__(seqs)

        if not self._validate_lengths():
            raise AlignmentError("All sequences need to be of equal length.")

        if score is not None:
            self._score = float(score)
        self._start_end_positions = start_end_positions

    def distances(self, distance_fn=None):
        """Compute distances between all pairs of sequences

        Parameters
        ----------
        distance_fn : function, optional
            Function for computing the distance between a pair of sequences.
            This must take two sequences as input (as `skbio.Sequence` objects)
            and return a single integer or float value. Defaults to the default
            distance function used by `skbio.Sequence.distance`.

        Returns
        -------
        skbio.DistanceMatrix
            Matrix containing the distances between all pairs of sequences.

        See Also
        --------
        skbio.Sequence.distance

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
        >>> seqs = [DNA("A-CCGGG", id="s1"),
        ...         DNA("ATCC--G", id="s2"),
        ...         DNA("ATCCGGA", id="s3")]
        >>> a1 = Alignment(seqs)
        >>> print(a1.distances())
        3x3 distance matrix
        IDs:
        's1', 's2', 's3'
        Data:
        [[ 0.          0.42857143  0.28571429]
         [ 0.42857143  0.          0.42857143]
         [ 0.28571429  0.42857143  0.        ]]

        """
        return super(Alignment, self).distances(distance_fn)

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
            The list of sequence start/end positions, or ``None`` if this was
            not provided on object construction.

        Notes
        -----
        The start/end positions indicate the range of the unaligned sequences
        in the alignment. For example, if local alignment were performed on the
        sequences ACA and TACAT, depending on the specific algorithm that was
        used to perform the alignment, the start/end positions would likely be:
        ``[(0,2), (1,3)]``. This indicates that the first and last positions of
        the second sequence were not included in the alignment, and the
        aligned sequences were therefore:
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
            A list of position indices to be retained in the resulting
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
        >>> from skbio import Alignment
        >>> from skbio import DNA
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
        # indices to keep
        indices = [
            i for i in range(self.sequence_length()) if keep_position(i)]
        # iterate over sequences
        for sequence_index, seq in enumerate(self):
            # determine if we're keeping the current sequence
            if keep_seq(sequence_index, seq.id):
                # slice the current sequence with the indices
                result.append(seq[indices])
            # if we're not keeping the current sequence, move on to the next
            else:
                continue
        # pack the result up in the same type of object as the current object
        # and return it
        return self.__class__(result)

    def iter_positions(self, constructor=None):
        """Generator of Alignment positions (i.e., columns)

        Parameters
        ----------
        constructor : type, optional
            Constructor function for creating the positional values. By
            default, these will be the same type as corresponding
            `skbio.Sequence` in the `Alignment` object, but
            you can pass a `skbio.Sequence` class here to ensure that they are
            all of consistent type, or ``str`` to have them returned as
            strings.

        Returns
        -------
        GeneratorType
            Generator of lists of positional values in the `Alignment`
            (effectively the transpose of the alignment).

        See Also
        --------
        iter

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
        >>> sequences = [DNA('ACCGT--', id="seq1"),
        ...              DNA('AACCGGT', id="seq2")]
        >>> a1 = Alignment(sequences)
        >>> for position in a1.iter_positions():
        ...     print(position)
        [DNA('A', length=1, id='seq1'), DNA('A', length=1, id='seq2')]
        [DNA('C', length=1, id='seq1'), DNA('A', length=1, id='seq2')]
        [DNA('C', length=1, id='seq1'), DNA('C', length=1, id='seq2')]
        [DNA('G', length=1, id='seq1'), DNA('C', length=1, id='seq2')]
        [DNA('T', length=1, id='seq1'), DNA('G', length=1, id='seq2')]
        [DNA('-', length=1, id='seq1'), DNA('G', length=1, id='seq2')]
        [DNA('-', length=1, id='seq1'), DNA('T', length=1, id='seq2')]

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

    def majority_consensus(self):
        """Return the majority consensus sequence for the alignment.

        Returns
        -------
        skbio.Sequence
            The consensus sequence of the `Alignment`. In other words, at each
            position the most common character is chosen, and those characters
            are combined to create a new sequence. The sequence will not have
            its ID, description, or quality set; only the sequence will be set.
            The type of biological sequence that is returned will be the same
            type as the first sequence in the alignment, or ``Sequence`` if the
            alignment is empty.

        Notes
        -----
        If there are two characters that are equally abundant in the sequence
        at a given position, the choice of which of those characters will be
        present at that position in the result is arbitrary.

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> a1.majority_consensus()
        DNA('AT-C', length=4)

        """
        if self.is_empty():
            seq_constructor = Sequence
        else:
            seq_constructor = self[0].__class__

        # Counter.most_common returns an ordered list of the n most common
        # (sequence, count) items in Counter. Here we set n=1, and take only
        # the character, not the count.
        return seq_constructor(''.join(c.most_common(1)[0][0]
                               for c in self.position_counters()))

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
            than (or equal to) `maximum_gap_frequency` fraction of the
            sequences.

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
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
        gap_chars = self[0].gap_chars

        positions_to_keep = []
        for i, f in enumerate(position_frequencies):
            gap_frequency = sum([f[c] for c in gap_chars])
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
            than (or equal to) `maximum_gap_frequency` fraction of the
            positions.

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
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

        base_frequencies = self.kmer_frequencies(k=1, relative=True)
        gap_chars = self[0].gap_chars
        seqs_to_keep = []
        for seq, f in zip(self, base_frequencies):
            gap_frequency = sum([f[c] for c in gap_chars])
            if gap_frequency <= maximum_gap_frequency:
                seqs_to_keep.append(seq.id)
        return self.subalignment(seqs_to_keep=seqs_to_keep)

    def position_counters(self):
        """Return counts of characters at each position in the alignment

        Returns
        -------
        list
            List of ``collections.Counter`` objects, one for each position in
            the `Alignment`.

        See Also
        --------
        position_frequencies
        position_entropies

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
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
        kmer_frequencies

        Examples
        --------
        >>> from skbio import Alignment
        >>> from skbio import DNA
        >>> sequences = [DNA('AC--', id="seq1"),
        ...              DNA('AT-C', id="seq2"),
        ...              DNA('TT-C', id="seq3")]
        >>> a1 = Alignment(sequences)
        >>> position_freqs = a1.position_frequencies()
        >>> round(position_freqs[0]['A'], 3)
        0.667
        >>> round(position_freqs[1]['A'], 3)
        0.0

        """
        seq_count = self.sequence_count()
        result = []
        for pos_counter in self.position_counters():
            freqs = defaultdict(float)
            for char, count in viewitems(pos_counter):
                freqs[char] = count / seq_count
            result.append(freqs)
        return result

    def position_entropies(self, base=None,
                           nan_on_non_standard_chars=True):
        """Return Shannon entropy of positions in Alignment

        Parameters
        ----------
        base : float, optional
            Log base for entropy calculation. If not passed, default will be e
            (i.e., natural log will be computed).
        nan_on_non_standard_chars : bool, optional
            If True, the entropy at positions containing characters outside of
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
        >>> from skbio import Alignment
        >>> from skbio import DNA
        >>> sequences = [DNA('AA--', id="seq1"),
        ...              DNA('AC-C', id="seq2"),
        ...              DNA('AT-C', id="seq3"),
        ...              DNA('TG-C', id="seq4")]
        >>> a1 = Alignment(sequences)
        >>> print(a1.position_entropies())
        [0.56233514461880829, 1.3862943611198906, nan, nan]

        """
        result = []
        # handle empty Alignment case
        if self.is_empty():
            return result

        iupac_standard_characters = self[0].nondegenerate_chars
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
        >>> from skbio import Alignment
        >>> from skbio import DNA
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

    def sequence_logo(self, reverse=False):
        """Plots alignment as a sequence logo.

        Parameters
        ----------
        reverse: bool, optional
            If ``True``, letters on plot are ordered with the largest letters
            on the bottom, and the smallest on the top.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the sequence logo of the passed alignment.

        Examples
        --------
        ..plot::

            Define an alignment of DNA sequences with original ids (the names
            of the ids will not change the sequence logo, however Alignment
            requires original ids):

            >>> from skbio import Alignment, DNA
            >>> a1 = Alignment([DNA('ATT-GG-G', id='seq1'),
            ...                 DNA('ACGTT--G', id='seq2'),
            ...                 DNA('CTCTG--G', id='seq3'),
            ...                 DNA('GCT--T-G', id='seq4')])

            Plot the alignment as a sequence logo:

            >>> fig = a1.sequence_logo()

        """
        # catches an empty alignment so the user won't get an empty graph
        if self.sequence_count() == 0:
            raise AlignmentError("You have passed an empty Alignment.")

        # acquires amount of each nucleotide at a given point
        pos_freq = self.position_frequencies()

        # sets up graph with x and y ticks and labels
        fig, ax = plt.subplots()
        fig.set_size_inches([(len(pos_freq)*1.3), 6])
        ax.set_xticks(np.arange(0.35, len(pos_freq)))
        ax.set_xlim(right=len(pos_freq))
        ax.set_xticklabels(range(1, self.sequence_length()+1))
        ax.set_ylim(top=1.0, bottom=0.0)
        ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.0', '', '', '', '1.0'])

        # runs through every sequence and graphs it
        for idx, pos in enumerate(pos_freq):
            # checks for accepted gap characters
            for gap in BiologicalSequence.gap_alphabet():
                if gap in pos:
                    # instead of drawing an empty box for gaps, we delete gap
                    # positions
                    del pos[gap]
            self._draw_sequence_logo_position(ax, idx, pos, reverse)

        return fig

    def _draw_sequence_logo_position(self, ax, position, freqs, reverse):
        sorted_freqs = sorted(freqs.items(), key=operator.itemgetter(1),
                              reverse=reverse)
        dx = 0.7
        x = position
        y = 0
        for letter, dy in sorted_freqs:
            try:
                add_letter(ax, letter, x, y, dx, dy)
            except KeyError as e:
                raise AlignmentError("Character %r cannot currently be drawn "
                                     "in ``Alignment.sequence_logo``. Support "
                                     "for more characters will be added in"
                                     "the future." % e.args[0])
            y += dy

    def _validate_lengths(self):
        """Return ``True`` if all sequences same length, ``False`` otherwise
        """
        seq1_length = self.sequence_length()
        for seq in self:
            if seq1_length != len(seq):
                return False
        return True


class StockholmAlignment(Alignment):
    """Contains the metadata information in a Stockholm file alignment

    Parameters
    ----------
    seqs : list of `skbio.Sequence` objects
        The `skbio.Sequence` objects to load.
    gf : dict, optional
        GF info in the format {feature: info}
    gs : dict of dicts, optional
        GS info in the format {feature: {seqlabel: info}}
    gr : dict of dicts, optional
        GR info in the format {feature: {seqlabel: info}}
    gc : dict, optional
        GC info in the format {feature: info}

    Notes
    -----
    The Stockholm format is described in [1]_ and [2]_.

    If there are multiple references, include information for each R* line
    as a list, with reference 0 information in position 0 for all lists,
    etc. This list will be broken up into the appropriate bits for each
    reference on string formatting.

    If there are multiple trees included, use a list to store identifiers
    and trees, with position 0 holding identifier for tree in position 0,
    etc.

    References
    ----------
    .. [1] http://sonnhammer.sbc.su.se/Stockholm.html
    .. [2] http://en.wikipedia.org/wiki/Stockholm_format

    Examples
    --------
    Assume we have a basic stockholm file with the following contents::

        # STOCKHOLM 1.0
        seq1         ACC--G-GGGU
        seq2         UCC--G-GGGA
        #=GC SS_cons (((.....)))
        //

    >>> from skbio import RNA
    >>> from skbio.alignment import StockholmAlignment
    >>> from StringIO import StringIO
    >>> sto_in = StringIO("# STOCKHOLM 1.0\\n"
    ...                   "seq1     ACC--G-GGGU\\nseq2     UCC--G-GGGA\\n"
    ...                   "#=GC SS_cons (((.....)))\\n//")
    >>> sto_records = StockholmAlignment.from_file(sto_in, RNA)
    >>> sto = next(sto_records)
    >>> print(sto)
    # STOCKHOLM 1.0
    seq1          ACC--G-GGGU
    seq2          UCC--G-GGGA
    #=GC SS_cons  (((.....)))
    //
    >>> sto.gc
    {'SS_cons': '(((.....)))'}

    We can also write out information by instantiating the StockholmAlignment
    object and then printing it.

    >>> from skbio import RNA
    >>> from skbio.alignment import StockholmAlignment
    >>> seqs = [RNA("ACC--G-GGGU", id="seq1"),
    ...         RNA("UCC--G-GGGA", id="seq2")]
    >>> gf = {
    ... "RT": ["TITLE1",  "TITLE2"],
    ... "RA": ["Auth1;", "Auth2;"],
    ... "RL": ["J Mol Biol", "Cell"],
    ... "RM": ["11469857", "12007400"]}
    >>> sto = StockholmAlignment(seqs, gf=gf)
    >>> print(sto)
    # STOCKHOLM 1.0
    #=GF RN [1]
    #=GF RM 11469857
    #=GF RT TITLE1
    #=GF RA Auth1;
    #=GF RL J Mol Biol
    #=GF RN [2]
    #=GF RM 12007400
    #=GF RT TITLE2
    #=GF RA Auth2;
    #=GF RL Cell
    seq1          ACC--G-GGGU
    seq2          UCC--G-GGGA
    //
    """
    def __init__(self, seqs, gf=None, gs=None, gr=None, gc=None):
        self.gf = gf if gf else {}
        self.gs = gs if gs else {}
        self.gr = gr if gr else {}
        self.gc = gc if gc else {}
        super(StockholmAlignment, self).__init__(seqs)

    def __str__(self):
        """Parses StockholmAlignment into a string with stockholm format

        Returns
        -------
        str
            Stockholm formatted string containing all information in the object

        Notes
        -----
        If references are included in GF data, the RN lines are automatically
        generated if not provided.

        """

        # find length of leader info needed to make file pretty
        # 10 comes from the characters for '#=GF ' and the feature after label
        infolen = max(len(seq.id) for seq in self._data) + 10

        GF_lines = []
        GS_lines = []
        GC_lines = []
        # NOTE: EVERYTHING MUST BE COERECED TO STR in case int or float passed
        # add GF information if applicable
        if self.gf:
            skipfeatures = set(("NH", "RC", "RM", "RN", "RA", "RL"))
            for feature, value in self.gf.items():
                # list of features to skip and parse special later
                if feature in skipfeatures:
                    continue
                # list of features to parse special
                elif feature == "TN":
                    # trees must be in proper order of identifier then tree
                    ident = value if isinstance(value, list) else [value]
                    tree = self.gf["NH"] if isinstance(self.gf["NH"], list) \
                        else [self.gf["NH"]]
                    for ident, tree in zip(self.gf["TN"], self.gf["NH"]):
                        GF_lines.append(' '.join(["#=GF", "TN", str(ident)]))
                        GF_lines.append(' '.join(["#=GF", "NH", str(tree)]))
                elif feature == "RT":
                    # make sure each reference block stays together
                    # set up lists to zip in case some bits are missing
                    # create rn list if needed
                    default_none = [0]*len(value)
                    rn = self.gf.get("RN", ["[%i]" % x for x in
                                     range(1, len(value)+1)])
                    rm = self.gf.get("RM", default_none)
                    rt = self.gf.get("RT", default_none)
                    ra = self.gf.get("RA", default_none)
                    rl = self.gf.get("RL", default_none)
                    rc = self.gf.get("RC", default_none)
                    # order: RN, RM, RT, RA, RL, RC
                    for n, m, t, a, l, c in zip(rn, rm, rt, ra, rl, rc):
                        GF_lines.append(' '.join(["#=GF", "RN", n]))
                        if m:
                            GF_lines.append(' '.join(["#=GF", "RM", str(m)]))
                        if t:
                            GF_lines.append(' '.join(["#=GF", "RT", str(t)]))
                        if a:
                            GF_lines.append(' '.join(["#=GF", "RA", str(a)]))
                        if l:
                            GF_lines.append(' '.join(["#=GF", "RL", str(l)]))
                        if c:
                            GF_lines.append(' '.join(["#=GF", "RC", str(c)]))
                else:
                    # normal addition for everything else
                    if not isinstance(value, list):
                        value = [value]
                    for val in value:
                        GF_lines.append(' '.join(["#=GF", feature, str(val)]))

        # add GS information if applicable
        if self.gs:
            for feature in self.gs:
                for seqname in self.gs[feature]:
                    GS_lines.append(' '.join(["#=GS", seqname, feature,
                                             str(self.gs[feature][seqname])]))

        # add GC information if applicable
        if self.gc:
            for feature, value in viewitems(self.gc):
                leaderinfo = ' '.join(["#=GC", feature])
                spacer = ' ' * (infolen - len(leaderinfo))
                GC_lines.append(spacer.join([leaderinfo,
                                             str(self.gc[feature])]))

        sto_lines = ["# STOCKHOLM 1.0"] + GF_lines + GS_lines
        # create seq output along with GR info if applicable
        for label, seq in self.iteritems():
            spacer = ' ' * (infolen - len(label))
            sto_lines.append(spacer.join([label, str(seq)]))
            # GR info added for sequence
            for feature in viewkeys(self.gr):
                value = self.gr[feature][label]
                leaderinfo = ' '.join(['#=GR', label, feature])
                spacer = ' ' * (infolen - len(leaderinfo))
                sto_lines.append(spacer.join([leaderinfo, value]))

        sto_lines.extend(GC_lines)
        # add final slashes to end of file
        sto_lines.append('//')

        return '\n'.join(sto_lines)

    def to_file(self, out_f):
        r"""Save the alignment to file in text format.

        Parameters
        ----------
        out_f : file-like object or filename
            File-like object to write serialized data to, or name of
            file. If it's a file-like object, it must have a ``write``
            method, and it won't be closed. Else, it is opened and
            closed after writing.

        See Also
        --------
        from_file
        """
        with open_file(out_f, 'w') as out_f:
            out_f.write(self.__str__())

    @staticmethod
    def _parse_gf_info(lines):
        """Takes care of parsing GF lines in stockholm plus special cases"""
        parsed = defaultdict(list)
        # needed for making each multi-line RT and NH one string
        rt = []
        nh = []
        lastline = ""
        for line in lines:
            try:
                init, feature, content = line.split(None, 2)
            except ValueError:
                raise StockholmParseError("Malformed GF line encountered!"
                                          "\n%s" % line.split(None, 2))
            if init != "#=GF":
                raise StockholmParseError("Non-GF line encountered!")

            # take care of adding multiline RT to the parsed information
            if lastline == "RT" and feature != "RT":
                # add rt line to the parsed dictionary
                rtline = " ".join(rt)
                rt = []
                parsed["RT"].append(rtline)
            elif feature == "RT":
                rt.append(content)
                lastline = feature
                continue

            # Take care of adding multiline NH to the parsed dictionary
            elif lastline == "NH" and feature != "NH":
                nhline = " ".join(nh)
                nh = []
                parsed["NH"].append(nhline)
            elif feature == "NH":
                nh.append(content)
                lastline = feature
                continue

            # add current feature to the parsed information
            parsed[feature].append(content)
            lastline = feature

        # removing unneccessary lists from parsed. Use .items() for py3 support
        for feature, value in parsed.items():
            # list of multi-line features to join into single string if needed
            if feature in ["CC"]:
                parsed[feature] = ' '.join(value)
            elif len(parsed[feature]) == 1:
                parsed[feature] = value[0]
        return parsed

    @staticmethod
    def _parse_gc_info(lines, strict=False, seqlen=-1):
        """Takes care of parsing GC lines in stockholm format"""
        parsed = {}
        for line in lines:
            try:
                init, feature, content = line.split(None, 2)
            except ValueError:
                raise StockholmParseError("Malformed GC line encountered!\n%s"
                                          % line.split(None, 2))
            if init != "#=GC":
                raise StockholmParseError("Non-GC line encountered!")

            # add current feature to the parsed information
            if feature in parsed:
                if strict:
                    raise StockholmParseError("Should not have multiple lines "
                                              "with the same feature: %s" %
                                              feature)
            else:
                parsed[feature] = [content]

        # removing unneccessary lists from parsed. Use .items() for py3 support
        for feature, value in parsed.items():
            parsed[feature] = ''.join(value)
            if strict:
                if len(value) != seqlen:
                    raise StockholmParseError("GC must have exactly one char "
                                              "per position in alignment!")

        return parsed

    @staticmethod
    def _parse_gs_gr_info(lines, strict=False, seqlen=-1):
        """Takes care of parsing GS and GR lines in stockholm format"""
        parsed = {}
        parsetype = ""
        for line in lines:
            try:
                init, label, feature, content = line.split(None, 3)
            except ValueError:
                raise StockholmParseError("Malformed GS/GR line encountered!"
                                          "\n%s" % line.split(None, 3))
            if parsetype == "":
                parsetype = init
            elif init != parsetype:
                    raise StockholmParseError("Non-GS/GR line encountered!")

            # parse each line, taking into account interleaved format
            if feature in parsed and label in parsed[feature]:
                # interleaved format, so need list of content
                parsed[feature][label].append(content)
            else:
                parsed[feature] = {label: [content]}

        # join all the crazy lists created during parsing
        for feature in parsed:
            for label, content in parsed[feature].items():
                parsed[feature][label] = ''.join(content)
                if strict:
                    if len(parsed[feature][label]) != seqlen:
                        raise StockholmParseError("GR must have exactly one "
                                                  "char per position in the "
                                                  "alignment!")
        return parsed

    @classmethod
    def from_file(cls, infile, seq_constructor, strict=False):
        r"""yields StockholmAlignment objects from a stockholm file.

        Parameters
        ----------
        infile : open file object
            An open stockholm file.

        seq_constructor : Sequence object
            The Sequence object that corresponds to what the
            stockholm file holds. See skbio

        strict : bool (optional)
            Turns on strict parsing of GR and GC lines to ensure one char per
             position. Default: False

        Returns
        -------
        Iterator of StockholmAlignment objects

        Raises
        ------
        skbio.StockholmParseError
            If any lines are found that don't conform to stockholm format
        """
        # make sure first line is corect
        line = infile.readline()
        if not line.startswith("# STOCKHOLM 1.0"):
            raise StockholmParseError("Incorrect header found")
        gs_lines = []
        gf_lines = []
        gr_lines = []
        gc_lines = []
        # OrderedDict used so sequences maintain same order as in file
        seqs = OrderedDict()
        for line in infile:
            line = line.strip()
            if line == "" or line.startswith("# S"):
                # skip blank lines or secondary headers
                continue
            elif line == "//":
                # parse the record since we are at its end
                # build the seuence list for alignment construction
                seqs = [seq_constructor(seq, id=_id) for _id, seq in
                        viewitems(seqs)]
                # get length of sequences in the alignment
                seqlen = len(seqs[0][1])

                # parse information lines
                gf = cls._parse_gf_info(gf_lines)
                gs = cls._parse_gs_gr_info(gs_lines)
                gr = cls._parse_gs_gr_info(gr_lines, strict, seqlen)
                gc = cls._parse_gc_info(gc_lines, strict, seqlen)

                # yield the actual stockholm object
                yield cls(seqs, gf, gs, gr, gc)

                # reset all storage variables
                gs_lines = []
                gf_lines = []
                gr_lines = []
                gc_lines = []
                seqs = OrderedDict()
            # add the metadata lines to the proper lists
            elif line.startswith("#=GF"):
                gf_lines.append(line)
            elif line.startswith("#=GS"):
                gs_lines.append(line)
            elif line.startswith("#=GR"):
                gr_lines.append(line)
            elif line.startswith("#=GC"):
                gc_lines.append(line)
            else:
                lineinfo = line.split()
                # assume sequence since nothing else in format is left
                # in case of interleaved format, need to do check
                if lineinfo[0] in seqs:
                    sequence = seqs[lineinfo[0]]
                    seqs[lineinfo[0]] = ''.join([sequence, lineinfo[1]])
                else:
                    seqs[lineinfo[0]] = lineinfo[1]
