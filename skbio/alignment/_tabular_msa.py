# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import collections
import copy

from future.builtins import range
from future.utils import viewkeys, viewvalues
import numpy as np
import pandas as pd

from skbio._base import SkbioObject, MetadataMixin, PositionalMetadataMixin
from skbio.sequence import Sequence
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util._decorator import experimental, classonlymethod, overrides
from skbio.util._misc import resolve_key

from skbio.alignment._repr import _TabularMSAReprBuilder


_Shape = collections.namedtuple('Shape', ['sequence', 'position'])


class TabularMSA(MetadataMixin, PositionalMetadataMixin, SkbioObject):
    """Store a multiple sequence alignment in tabular (row/column) form.

    Parameters
    ----------
    sequences : iterable of alphabet-aware skbio sequence objects, TabularMSA
        Aligned sequences in the MSA. Sequences must be the same type, length,
        and have an alphabet. For example, `sequences` could be an iterable of
        ``DNA``, ``RNA``, or ``Protein`` objects. If `sequences` is a
        ``TabularMSA``, its `metadata`, `positional_metadata`, and `index` will
        be used unless overridden by parameters `metadata`,
        `positional_metadata`, and `minter`/`index`.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire MSA. A shallow copy of
        the ``dict`` will be made.
    positional_metadata : pd.DataFrame consumable, optional
        Arbitrary metadata which applies to each position in the MSA. Must be
        able to be passed directly to ``pd.DataFrame`` constructor. Each column
        of metadata must be the same length as the number of positions in the
        MSA. A shallow copy of the positional metadata will be made.
    minter : callable or metadata key, optional
        If provided, defines an index label for each sequence in `sequences`.
        Can either be a callable accepting a single argument (each sequence) or
        a key into each sequence's ``metadata`` attribute.
    index : pd.Index consumable, optional
        Index containing labels for `sequences`. Must be the same length as
        `sequences`. Must be able to be passed directly to ``pd.Index``
        constructor.

    Raises
    ------
    ValueError
        If `minter` and `index` are both provided.
    ValueError
        If `index` is not the same length as `sequences`.

    See Also
    --------
    skbio.sequence.DNA
    skbio.sequence.RNA
    skbio.sequence.Protein
    pandas.DataFrame
    pandas.Index

    Notes
    -----
    If `minter` or `index` are not provided, default pandas labels will be
    used: integer labels ``0..(N-1)``, where ``N`` is the number of sequences.

    """

    @property
    @experimental(as_of='0.4.0-dev')
    def dtype(self):
        """Data type of the stored sequences.

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa.dtype
        <class 'skbio.sequence._dna.DNA'>
        >>> msa.dtype is DNA
        True

        """
        return type(self._get_sequence(0)) if len(self) > 0 else None

    @property
    @experimental(as_of='0.4.0-dev')
    def shape(self):
        """Number of sequences (rows) and positions (columns).

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA

        Create a ``TabularMSA`` object with 2 sequences and 3 positions:

        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa.shape
        Shape(sequence=2, position=3)
        >>> msa.shape == (2, 3)
        True

        Dimensions can be accessed by index or by name:

        >>> msa.shape[0]
        2
        >>> msa.shape.sequence
        2
        >>> msa.shape[1]
        3
        >>> msa.shape.position
        3

        """
        sequence_count = len(self)

        if sequence_count > 0:
            position_count = len(self._get_sequence(0))
        else:
            position_count = 0

        return _Shape(sequence=sequence_count, position=position_count)

    @property
    @experimental(as_of='0.4.0-dev')
    def index(self):
        """Index containing labels along the sequence axis.

        Returns
        -------
        pd.Index
            Index containing sequence labels.

        See Also
        --------
        reassign_index

        Notes
        -----
        This property can be set and deleted. Deleting the index will reset the
        index to the ``TabularMSA`` constructor's default.

        Examples
        --------
        Create a ``TabularMSA`` object with sequences labeled by sequence
        identifier:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'a'}),
        ...         DNA('AC-', metadata={'id': 'b'})]
        >>> msa = TabularMSA(seqs, minter='id')

        Retrieve index:

        >>> msa.index
        Index(['a', 'b'], dtype='object')

        Set index:

        >>> msa.index = ['seq1', 'seq2']
        >>> msa.index
        Index(['seq1', 'seq2'], dtype='object')

        Delete index:

        >>> del msa.index
        >>> msa.index
        Int64Index([0, 1], dtype='int64')

        """
        return self._seqs.index

    @index.setter
    def index(self, index):
        self._seqs.index = index

    @index.deleter
    def index(self):
        self.reassign_index()

    @classonlymethod
    @experimental(as_of="0.4.0-dev")
    def from_dict(cls, dictionary):
        """Create a ``TabularMSA`` from a ``dict``.

        Parameters
        ----------
        dictionary : dict
            Dictionary mapping keys to alphabet-aware scikit-bio sequence
            objects. The ``TabularMSA`` object will have its index labels set
            to the keys in the dictionary.

        Returns
        -------
        TabularMSA
            ``TabularMSA`` object constructed from the keys and sequences in
            `dictionary`.

        See Also
        --------
        to_dict
        sort

        Notes
        -----
        The order of sequences and index labels in the resulting ``TabularMSA``
        object is arbitrary. Use ``TabularMSA.sort`` to set a different order.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> seqs = {'a': DNA('ACGT'), 'b': DNA('A--T')}
        >>> msa = TabularMSA.from_dict(seqs)

        """
        # Python 2 and 3 guarantee same order of iteration as long as no
        # modifications are made to the dictionary between calls:
        #     https://docs.python.org/2/library/stdtypes.html#dict.items
        #     https://docs.python.org/3/library/stdtypes.html#
        #         dictionary-view-objects
        return cls(viewvalues(dictionary), index=viewkeys(dictionary))

    @experimental(as_of='0.4.0-dev')
    def __init__(self, sequences, metadata=None, positional_metadata=None,
                 minter=None, index=None):
        if isinstance(sequences, TabularMSA):
            if metadata is None and sequences.has_metadata():
                metadata = sequences.metadata
            if (positional_metadata is None and
                    sequences.has_positional_metadata()):
                positional_metadata = sequences.positional_metadata
            if minter is None and index is None:
                index = sequences.index

        self._seqs = pd.Series([])
        self.extend(sequences, minter=minter, index=index)

        MetadataMixin._init_(self, metadata=metadata)
        PositionalMetadataMixin._init_(
            self, positional_metadata=positional_metadata)

    @experimental(as_of='0.4.0-dev')
    def __repr__(self):
        pep8_line_length_limit = 79
        length_taken_by_docstring_indent = 8
        width = pep8_line_length_limit - length_taken_by_docstring_indent
        return _TabularMSAReprBuilder(
            obj=self,
            width=width,
            indent=4).build()

    def _repr_stats(self):
        return [("sequence count", str(self.shape.sequence)),
                ("position count", str(self.shape.position))]

    @experimental(as_of='0.4.0-dev')
    def __bool__(self):
        """Boolean indicating whether the MSA is empty or not.

        Returns
        -------
        bool
            ``False`` if there are no sequences, OR if there are no positions
            (i.e., all sequences are empty). ``True`` otherwise.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA

        MSA with sequences and positions:

        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> bool(msa)
        True

        No sequences:

        >>> msa = TabularMSA([])
        >>> bool(msa)
        False

        No positions:

        >>> msa = TabularMSA([DNA(''), DNA('')])
        >>> bool(msa)
        False

        """
        # It is impossible to have 0 sequences and >0 positions.
        return self.shape.position > 0

    # Python 2 compatibility.
    __nonzero__ = __bool__

    @experimental(as_of='0.4.0-dev')
    def __contains__(self, label):
        """Determine if an index label is in this MSA.

        Parameters
        ----------
        label : hashable
            Label to search for in this MSA.

        Returns
        -------
        bool
            Indicates whether `label` is in this MSA.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')], index=['l1', 'l2'])
        >>> 'l1' in msa
        True
        >>> 'l2' in msa
        True
        >>> 'l3' in msa
        False

        """
        return label in self.index

    @experimental(as_of='0.4.0-dev')
    def __len__(self):
        """Number of sequences in the MSA.

        Returns
        -------
        int
            Number of sequences in the MSA (i.e., size of the 1st dimension).

        Notes
        -----
        This is equivalent to ``msa.shape[0]``.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> len(msa)
        2
        >>> msa = TabularMSA([])
        >>> len(msa)
        0
        """
        return len(self._seqs)

    @experimental(as_of='0.4.0-dev')
    def __iter__(self):
        """Iterate over sequences in the MSA.

        Yields
        ------
        alphabet-aware scikit-bio sequence object
            Each sequence in the order they are stored in the MSA.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> for seq in msa:
        ...     str(seq)
        'ACG'
        'AC-'

        """
        return iter(self._seqs)

    @experimental(as_of='0.4.0-dev')
    def __reversed__(self):
        """Iterate in reverse order over sequences in the MSA.

        Yields
        ------
        alphabet-aware scikit-bio sequence object
            Each sequence in reverse order from how they are stored in the MSA.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> for seq in reversed(msa):
        ...     str(seq)
        'AC-'
        'ACG'

        """
        return reversed(self._seqs)

    @experimental(as_of='0.4.0-dev')
    def __str__(self):
        return self.__repr__()

    @experimental(as_of='0.4.0-dev')
    def __eq__(self, other):
        """Determine if this MSA is equal to another.

        ``TabularMSA`` objects are equal if their sequences, index, metadata,
        and positional metadata are equal.

        Parameters
        ----------
        other : TabularMSA
            MSA to test for equality against.

        Returns
        -------
        bool
            Indicates whether this MSA is equal to `other`.

        Examples
        --------
        >>> from skbio import DNA, RNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa == msa
        True

        MSAs with different sequence characters are not equal:

        >>> msa == TabularMSA([DNA('ACG'), DNA('--G')])
        False

        MSAs with different types of sequences (different ``dtype``) are not
        equal:

        >>> msa == TabularMSA([RNA('ACG'), RNA('AC-')])
        False

        MSAs with different sequence metadata are not equal:

        >>> msa == TabularMSA([DNA('ACG', metadata={'id': 'a'}), DNA('AC-')])
        False

        MSAs with different index labels are not equal:

        >>> msa == TabularMSA([DNA('ACG'), DNA('AC-')], minter=str)
        False

        MSAs with different metadata are not equal:

        >>> msa == TabularMSA([DNA('ACG'), DNA('AC-')],
        ...                   metadata={'id': 'msa-id'})
        False

        MSAs with different positional metadata are not equal:

        >>> msa == TabularMSA([DNA('ACG'), DNA('AC-')],
        ...                   positional_metadata={'prob': [3, 2, 1]})
        False

        """
        if not isinstance(other, TabularMSA):
            return False

        if not MetadataMixin._eq_(self, other):
            return False

        if not PositionalMetadataMixin._eq_(self, other):
            return False

        return self._seqs.equals(other._seqs)

    @experimental(as_of='0.4.0-dev')
    def __ne__(self, other):
        """Determine if this MSA is not equal to another.

        ``TabularMSA`` objects are not equal if their sequences, index,
        metadata, or positional metadata are not equal.

        Parameters
        ----------
        other : TabularMSA
            MSA to test for inequality against.

        Returns
        -------
        bool
            Indicates whether this MSA is not equal to `other`.

        See Also
        --------
        __eq__

        """
        return not (self == other)

    @experimental(as_of='0.4.0-dev')
    def __copy__(self):
        """Return a shallow copy of this MSA.

        Returns
        -------
        TabularMSA
            Shallow copy of this MSA. Sequence objects will be shallow-copied.

        See Also
        --------
        __deepcopy__

        """
        seqs = (copy.copy(seq) for seq in self._seqs)

        # Copying index isn't necessary because pd.Index is immutable.
        msa_copy = self.__class__(sequences=seqs, index=self.index,
                                  metadata=None,
                                  positional_metadata=None)

        msa_copy._metadata = MetadataMixin._copy_(self)
        msa_copy._positional_metadata = PositionalMetadataMixin._copy_(self)

        return msa_copy

    @experimental(as_of='0.4.0-dev')
    def __deepcopy__(self, memo):
        """Return a deep copy of this MSA.

        Returns
        -------
        TabularMSA
            Deep copy of this MSA.

        See Also
        --------
        __copy__

        """
        seqs = (copy.deepcopy(seq, memo) for seq in self._seqs)

        # Copying index isn't necessary because pd.Index is immutable.
        msa_copy = self.__class__(sequences=seqs, index=self.index,
                                  metadata=None, positional_metadata=None)

        msa_copy._metadata = MetadataMixin._deepcopy_(self, memo)
        msa_copy._positional_metadata = \
            PositionalMetadataMixin._deepcopy_(self, memo)

        return msa_copy

    @experimental(as_of='0.4.0-dev')
    def iter_positions(self, reverse=False):
        """Iterate over positions (columns) in the MSA.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, iterate over positions in reverse order.

        Yields
        ------
        Sequence
            Each position in the order they are stored in the MSA.

        See Also
        --------
        __iter__
        __reversed__
        skbio.sequence.Sequence.concat

        Notes
        -----
        Each position will be yielded as *exactly* a ``Sequence`` object,
        regardless of this MSA's ``dtype``. ``Sequence`` is used because a
        position is an artifact of multiple sequence alignment and is not a
        real biological sequence.

        Each ``Sequence`` object will have its corresponding MSA positional
        metadata stored as ``metadata``.

        Sequences will have their positional metadata concatenated using an
        outer join. See ``Sequence.concat(how='outer')`` for details.

        Examples
        --------
        Create an MSA with positional metadata:

        >>> from skbio import DNA, TabularMSA
        >>> sequences = [DNA('ACG'),
        ...              DNA('A-T')]
        >>> msa = TabularMSA(sequences,
        ...                  positional_metadata={'prob': [3, 1, 2]})

        Iterate over positions:

        >>> for position in msa.iter_positions():
        ...     position
        ...     print()
        Sequence
        -------------
        Metadata:
            'prob': 3
        Stats:
            length: 2
        -------------
        0 AA
        <BLANKLINE>
        Sequence
        -------------
        Metadata:
            'prob': 1
        Stats:
            length: 2
        -------------
        0 C-
        <BLANKLINE>
        Sequence
        -------------
        Metadata:
            'prob': 2
        Stats:
            length: 2
        -------------
        0 GT
        <BLANKLINE>

        Note that MSA positional metadata is stored as ``metadata`` on each
        ``Sequence`` object.

        Iterate over positions in reverse order:

        >>> for position in msa.iter_positions(reverse=True):
        ...     position
        ...     print('')
        Sequence
        -------------
        Metadata:
            'prob': 2
        Stats:
            length: 2
        -------------
        0 GT
        <BLANKLINE>
        Sequence
        -------------
        Metadata:
            'prob': 1
        Stats:
            length: 2
        -------------
        0 C-
        <BLANKLINE>
        Sequence
        -------------
        Metadata:
            'prob': 3
        Stats:
            length: 2
        -------------
        0 AA
        <BLANKLINE>

        """
        indices = range(self.shape.position)
        if reverse:
            indices = reversed(indices)

        return (self._get_position(index) for index in indices)

    @experimental(as_of='0.4.0-dev')
    def consensus(self):
        """Compute the majority consensus sequence for this MSA.

        The majority consensus sequence contains the most common character at
        each position in this MSA. Ties will be broken in an arbitrary manner.

        Returns
        -------
        Sequence
            The majority consensus sequence for this MSA. The type of sequence
            returned will be the same as this MSA's ``dtype`` or ``Sequence``
            if this MSA does not contain any sequences. The majority consensus
            sequence will have its positional metadata set to this MSA's
            positional metadata if present.

        Notes
        -----
        The majority consensus sequence will use this MSA's default gap
        character (``dtype.default_gap_char``) to represent gap majority at a
        position, regardless of the gap characters present at that position.

        Different gap characters at a position are **not** treated as distinct
        characters. All gap characters at a position contribute to that
        position's gap consensus.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> sequences = [DNA('AC---'),
        ...              DNA('AT-C.'),
        ...              DNA('TT-CG')]
        >>> msa = TabularMSA(sequences,
        ...                  positional_metadata={'prob': [2, 1, 2, 3, 5]})
        >>> msa.consensus()
        DNA
        -----------------------------
        Positional metadata:
            'prob': <dtype: int64>
        Stats:
            length: 5
            has gaps: True
            has degenerates: False
            has non-degenerates: True
            GC-content: 33.33%
        -----------------------------
        0 AT-C-

        Note that the last position in the MSA has more than one type of gap
        character. These are not treated as distinct characters; both types of
        gap characters contribute to the position's consensus. Also note that
        ``DNA.default_gap_char`` is used to represent gap majority at a
        position (``'-'``).

        """
        dtype = self.dtype
        if dtype is None:
            dtype = Sequence

        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata

        consensus = []
        for position in self.iter_positions():
            freqs = position.frequencies()

            gap_freq = 0
            for gap_char in dtype.gap_chars:
                if gap_char in freqs:
                    gap_freq += freqs.pop(gap_char)
            assert dtype.default_gap_char not in freqs
            freqs[dtype.default_gap_char] = gap_freq

            consensus.append(collections.Counter(freqs).most_common(1)[0][0])

        return dtype(''.join(consensus),
                     positional_metadata=positional_metadata)

    @experimental(as_of='0.4.0-dev')
    def gap_frequencies(self, axis='sequence', relative=False):
        """Compute frequency of gap characters across an axis.

        Parameters
        ----------
        axis : {'sequence', 'position'}, optional
            Axis to compute gap character frequencies across. If 'sequence' or
            0, frequencies are computed for each position in the MSA. If
            'position' or 1, frequencies are computed for each sequence.
        relative : bool, optional
            If ``True``, return the relative frequency of gap characters
            instead of the count.

        Returns
        -------
        1D np.ndarray (int or float)
            Vector of gap character frequencies across the specified axis. Will
            have ``int`` dtype if ``relative=False`` and ``float`` dtype if
            ``relative=True``.

        Raises
        ------
        ValueError
            If `axis` is invalid.

        Notes
        -----
        If there are no positions in the MSA, ``axis='position'``, **and**
        ``relative=True``, the relative frequency of gap characters in each
        sequence will be ``np.nan``.

        Examples
        --------
        Compute frequency of gap characters for each position in the MSA (i.e.,
        *across* the sequence axis):

        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'),
        ...                   DNA('A--'),
        ...                   DNA('AC.'),
        ...                   DNA('AG.')])
        >>> msa.gap_frequencies()
        array([0, 1, 3])

        Compute relative frequencies across the same axis:

        >>> msa.gap_frequencies(relative=True)
        array([ 0.  ,  0.25,  0.75])

        Compute frequency of gap characters for each sequence (i.e., *across*
        the position axis):

        >>> msa.gap_frequencies(axis='position')
        array([0, 2, 1, 1])

        """
        if self._is_sequence_axis(axis):
            seq_iterator = self.iter_positions()
            length = self.shape.sequence
        else:
            seq_iterator = self
            length = self.shape.position

        gap_freqs = []
        for seq in seq_iterator:
            # Not using Sequence.frequencies(relative=relative) because each
            # gap character's relative frequency is computed separately and
            # must be summed. This is less precise than summing the absolute
            # frequencies of gap characters and dividing by the length. Likely
            # not a big deal for typical gap characters ('-', '.') but can be
            # problematic as the number of gap characters grows (we aren't
            # guaranteed to always have two gap characters). See unit tests for
            # an example.
            freqs = seq.frequencies(chars=self.dtype.gap_chars)
            gap_freqs.append(sum(viewvalues(freqs)))

        gap_freqs = np.asarray(gap_freqs, dtype=float if relative else int)

        if relative:
            gap_freqs /= length

        return gap_freqs

    @experimental(as_of='0.4.0-dev')
    def reassign_index(self, mapping=None, minter=None):
        """Reassign index labels to sequences in this MSA.

        Parameters
        ----------
        mapping : dict-like or callable, optional
            Dictionary or callable that maps existing labels to new labels. Any
            label without a mapping will remain the same.
        minter : callable or metadata key, optional
            If provided, defines an index label for each sequence. Can either
            be a callable accepting a single argument (each sequence) or a key
            into each sequence's ``metadata`` attribute.

        Raises
        ------
        ValueError
            If `mapping` and `minter` are both provided.

        See Also
        --------
        index

        Notes
        -----
        If neither `mapping` nor `minter` are provided, default pandas labels
        will be used: integer labels ``0..(N-1)``, where ``N`` is the number of
        sequences.

        Examples
        --------
        Create a ``TabularMSA`` object with default index labels:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'a'}),
        ...         DNA('AC-', metadata={'id': 'b'})]
        >>> msa = TabularMSA(seqs)
        >>> msa.index
        Int64Index([0, 1], dtype='int64')

        Assign new index to the MSA using each sequence's ID as a label:

        >>> msa.reassign_index(minter='id')
        >>> msa.index
        Index(['a', 'b'], dtype='object')

        Assign default index:

        >>> msa.reassign_index()
        >>> msa.index
        Int64Index([0, 1], dtype='int64')

        Alternatively, a mapping of existing labels to new labels may be passed
        via `mapping`:

        >>> msa.reassign_index(mapping={0: 'seq1', 1: 'seq2'})
        >>> msa.index
        Index(['seq1', 'seq2'], dtype='object')

        """
        if mapping is not None and minter is not None:
            raise ValueError(
                "Cannot use both `mapping` and `minter` at the same time.")
        if mapping is not None:
            self._seqs.rename(mapping, inplace=True)
        elif minter is not None:
            index = [resolve_key(seq, minter) for seq in self._seqs]

            # Cast to Index to identify tuples as a MultiIndex to match
            # pandas constructor. Just setting would make an index of tuples.
            self.index = pd.Index(index)
        else:
            self._seqs.reset_index(drop=True, inplace=True)

    @experimental(as_of='0.4.0-dev')
    def append(self, sequence, minter=None, index=None):
        """Append a sequence to the MSA without recomputing alignment.

        Parameters
        ----------
        sequence : alphabet-aware scikit-bio sequence object
            Sequence to be appended. Must match the dtype of the MSA and the
            number of positions in the MSA.
        minter : callable or metadata key, optional
            Used to create an index label for the sequence being appended. If
            callable, it generates a label directly. Otherwise it's treated as
            a key into the sequence metadata. Note that `minter` cannot be
            combined with `index`.
        index : object, optional
            Index label to use for the appended sequence. Note that `index`
            cannot be combined with `minter`.

        Raises
        ------
        ValueError
            If both `minter` and `index` are provided.
        ValueError
            If neither `minter` nor `index` are provided and the MSA has a
            non-default index.
        TypeError
            If the sequence object is a type that doesn't have an alphabet.
        TypeError
            If the type of the sequence does not match the dtype of the MSA.
        ValueError
            If the length of the sequence does not match the number of
            positions in the MSA.

        See Also
        --------
        extend
        reassign_index

        Notes
        -----
        If neither `minter` nor `index` are provided and this MSA has default
        index labels, the new index label will be auto-incremented.

        The MSA is not automatically re-aligned when a sequence is appended.
        Therefore, this operation is not necessarily meaningful on its own.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACGT')])
        >>> msa.append(DNA('AG-T'))
        >>> msa == TabularMSA([DNA('ACGT'), DNA('AG-T')])
        True

        Auto-incrementing index labels:

        >>> msa.index
        Int64Index([0, 1], dtype='int64')
        >>> msa.append(DNA('ACGA'))
        >>> msa.index
        Int64Index([0, 1, 2], dtype='int64')

        """
        if index is not None:
            index = [index]
        self.extend([sequence], minter=minter, index=index)

    @experimental(as_of='0.4.0-dev')
    def extend(self, sequences, minter=None, index=None):
        """Extend this MSA with sequences without recomputing alignment.

        Parameters
        ----------
        sequences : iterable of alphabet-aware scikit-bio sequence objects
            Sequences to be appended. Must match the dtype of the MSA and the
            number of positions in the MSA.
        minter : callable or metadata key, optional
            Used to create index labels for the sequences being appended. If
            callable, it generates a label directly. Otherwise it's treated as
            a key into the sequence metadata. Note that `minter` cannot be
            combined with `index`.
        index : pd.Index consumable, optional
            Index labels to use for the appended sequences. Must be the same
            length as `sequences`. Must be able to be passed directly to
            ``pd.Index`` constructor. Note that `index` cannot be combined
            with `minter`.

        Raises
        ------
        ValueError
            If both `minter` and `index` are both provided.
        ValueError
            If neither `minter` nor `index` are provided and the MSA has a
            non-default index.
        ValueError
            If `index` is not the same length as `sequences`.
        TypeError
            If `sequences` contains a type that does not have an alphabet.
        TypeError
            If `sequence` contains a type that does not match the dtype of the
            MSA.
        ValueError
            If the length of a sequence does not match the number of positions
            in the MSA.

        See Also
        --------
        append
        reassign_index

        Notes
        -----
        If neither `minter` nor `index` are provided and this MSA has default
        index labels, the new index labels will be auto-incremented.

        The MSA is not automatically re-aligned when appending sequences.
        Therefore, this operation is not necessarily meaningful on its own.

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACGT')])
        >>> msa.extend([DNA('AG-T'), DNA('-G-T')])
        >>> msa == TabularMSA([DNA('ACGT'), DNA('AG-T'), DNA('-G-T')])
        True

        Auto-incrementing index labels:

        >>> msa.index
        Int64Index([0, 1, 2], dtype='int64')
        >>> msa.extend([DNA('ACGA'), DNA('AC-T'), DNA('----')])
        >>> msa.index
        Int64Index([0, 1, 2, 3, 4, 5], dtype='int64')

        """
        if minter is not None and index is not None:
            raise ValueError(
                "Cannot use both `minter` and `index` at the same time.")

        sequences = list(sequences)

        if minter is None and index is None:
            if self.index.equals(pd.Index(np.arange(len(self)))):
                index = range(len(self), len(self) + len(sequences))
            else:
                raise ValueError(
                    "MSA does not have default index labels, must provide "
                    "a `minter` or `index` for sequence(s).")
        elif minter is not None:
            index = [resolve_key(seq, minter) for seq in sequences]

        # Cast to Index to identify tuples as a MultiIndex to match
        # pandas constructor. Just setting would make an index of tuples.
        if not isinstance(index, pd.Index):
            index = pd.Index(index)

        self._assert_valid_sequences(sequences)

        # pandas doesn't give a user-friendly error message if we pass through.
        if len(sequences) != len(index):
            raise ValueError(
                "Number of sequences (%d) must match index length (%d)" %
                (len(sequences), len(index)))
        self._seqs = self._seqs.append(pd.Series(sequences, index=index))

    def _assert_valid_sequences(self, sequences):
        if not sequences:
            return

        if len(self):
            expected_dtype = self.dtype
            expected_length = self.shape.position
        else:
            sequence = sequences[0]
            expected_dtype = type(sequence)
            if not issubclass(expected_dtype, IUPACSequence):
                raise TypeError(
                    "Each sequence must be a scikit-bio sequence object "
                    "that has an alphabet, not type %r" %
                    expected_dtype.__name__)
            expected_length = len(sequence)

        for sequence in sequences:
            dtype = type(sequence)
            if dtype is not expected_dtype:
                raise TypeError(
                    "Sequences in MSA must have matching type. Type %r does "
                    "not match type %r" % (dtype.__name__,
                                           expected_dtype.__name__))

            length = len(sequence)
            if length != expected_length:
                raise ValueError(
                    "Each sequence's length must match the number of "
                    "positions in the MSA: %d != %d"
                    % (length, expected_length))

    def join(self, other, how='strict'):
        """Join this MSA with another by sequence (horizontally).

        Sequences will be joined by index labels. MSA ``positional_metadata``
        will be joined by columns. Use `how` to control join behavior.

        Alignment is **not** recomputed during join operation (see *Notes*
        section for details).

        Parameters
        ----------
        other : TabularMSA
            MSA to join with. Must have same ``dtype`` as this MSA.
        how : {'strict', 'inner', 'outer', 'left', 'right'}, optional
            How to join the sequences and MSA `positional_metadata`:

            * ``'strict'``: MSA indexes and `positional_metadata` columns must
              match

            * ``'inner'``: an inner-join of the MSA indexes and
              ``positional_metadata`` columns (only the shared set of index
              labels and columns are used)

            * ``'outer'``: an outer-join of the MSA indexes and
              ``positional_metadata`` columns (all index labels and columns are
              used). Unshared sequences will be padded with the MSA's default
              gap character (``TabularMSA.dtype.default_gap_char``). Unshared
              columns will be padded with NaN.

            * ``'left'``: a left-outer-join of the MSA indexes and
              ``positional_metadata`` columns (this MSA's index labels and
              columns are used). Padding of unshared data is handled the same
              as ``'outer'``.

            * ``'right'``: a right-outer-join of the MSA indexes and
              ``positional_metadata`` columns (`other` index labels and columns
              are used). Padding of unshared data is handled the same as
              ``'outer'``.

        Returns
        -------
        TabularMSA
            Joined MSA. There is no guaranteed ordering to its index (call
            ``sort`` to define one).

        Raises
        ------
        ValueError
            If `how` is invalid.
        ValueError
            If either the index of this MSA or the index of `other` contains
            duplicates.
        ValueError
            If ``how='strict'`` and this MSA's index doesn't match with
            `other`.
        ValueError
            If ``how='strict'`` and this MSA's ``positional_metadata`` columns
            don't match with `other`.
        TypeError
            If `other` is not a subclass of ``TabularMSA``.
        TypeError
            If the ``dtype`` of `other` does not match this MSA's ``dtype``.

        See Also
        --------
        extend
        sort
        skbio.sequence.Sequence.concat

        Notes
        -----
        The join operation does not automatically perform re-alignment;
        sequences are simply joined together. Therefore, this operation is not
        necessarily meaningful on its own.

        The index labels of this MSA must be unique. Likewise, the index labels
        of `other` must be unique.

        The MSA-wide and per-sequence metadata (``TabularMSA.metadata`` and
        ``Sequence.metadata``) are not retained on the joined ``TabularMSA``.

        The positional metadata of the sequences will be outer-joined,
        regardless of `how` (using ``Sequence.concat(how='outer')``).

        If the join operation results in a ``TabularMSA`` without any
        sequences, the MSA's ``positional_metadata`` will not be set.

        Examples
        --------
        Join MSAs by sequence:

        >>> from skbio import DNA, TabularMSA
        >>> msa1 = TabularMSA([DNA('AC'),
        ...                    DNA('A-')])
        >>> msa2 = TabularMSA([DNA('G-T'),
        ...                    DNA('T--')])
        >>> joined = msa1.join(msa2)
        >>> joined == TabularMSA([DNA('ACG-T'),
        ...                       DNA('A-T--')])
        True

        Sequences are joined based on MSA index labels:

        >>> msa1 = TabularMSA([DNA('AC'),
        ...                    DNA('A-')], index=['a', 'b'])
        >>> msa2 = TabularMSA([DNA('G-T'),
        ...                    DNA('T--')], index=['b', 'a'])
        >>> joined = msa1.join(msa2)
        >>> joined == TabularMSA([DNA('ACT--'),
        ...                       DNA('A-G-T')], index=['a', 'b'])
        True

        By default both MSA indexes must match. Use ``how`` to specify an inner
        join:

        >>> msa1 = TabularMSA([DNA('AC'),
        ...                    DNA('A-'),
        ...                    DNA('-C')], index=['a', 'b', 'c'],
        ...                   positional_metadata={'col1': [42, 43],
        ...                                        'col2': [1, 2]})
        >>> msa2 = TabularMSA([DNA('G-T'),
        ...                    DNA('T--'),
        ...                    DNA('ACG')], index=['b', 'a', 'z'],
        ...                   positional_metadata={'col2': [3, 4, 5],
        ...                                        'col3': ['f', 'o', 'o']})
        >>> joined = msa1.join(msa2, how='inner')
        >>> joined == TabularMSA([DNA('A-G-T'),
        ...                       DNA('ACT--')], index=['b', 'a'],
        ...                      positional_metadata={'col2': [1, 2, 3, 4, 5]})
        True

        When performing an outer join (``'outer'``, ``'left'``, or
        ``'right'``), unshared sequences are padded with gaps and unshared
        ``positional_metadata`` columns are padded with NaN:

        >>> joined = msa1.join(msa2, how='outer')
        >>> joined == TabularMSA([DNA('ACT--'),
        ...                       DNA('A-G-T'),
        ...                       DNA('-C---'),
        ...                       DNA('--ACG')], index=['a', 'b', 'c', 'z'],
        ...                      positional_metadata={
        ...                          'col1': [42, 43, np.nan, np.nan, np.nan],
        ...                          'col2': [1, 2, 3, 4, 5],
        ...                          'col3': [np.nan, np.nan, 'f', 'o', 'o']})
        True

        """
        if how not in {'strict', 'inner', 'outer', 'left', 'right'}:
            raise ValueError(
                "`how` must be 'strict', 'inner', 'outer', 'left', or "
                "'right'.")

        self._assert_joinable(other)

        join_index, concat_kwargs = self._get_join_index(other, how)

        joined_seqs = []
        for label in join_index:
            left_seq = self._get_sequence_for_join(label)
            right_seq = other._get_sequence_for_join(label)

            joined_seqs.append(
                self.dtype.concat([left_seq, right_seq], how='outer'))

        # TODO: update when #1198 is implemented.
        joined_positional_metadata = None
        if joined_seqs:
            joined_positional_metadata = pd.concat(
                [self.positional_metadata, other.positional_metadata],
                ignore_index=True, **concat_kwargs)

            if not self.has_positional_metadata():
                del self.positional_metadata
            if not other.has_positional_metadata():
                del other.positional_metadata

        joined = self.__class__(joined_seqs, index=join_index,
                                positional_metadata=joined_positional_metadata)

        if not joined.has_positional_metadata():
            del joined.positional_metadata

        return joined

    def _assert_joinable(self, other):
        if not isinstance(other, TabularMSA):
            raise TypeError(
                "`other` must be a `TabularMSA` object, not type %r" %
                type(other).__name__)

        if self.dtype is not other.dtype:
            raise TypeError(
                "`other` dtype %r does not match this MSA's dtype %r" %
                (other.dtype if other.dtype is None else other.dtype.__name__,
                 self.dtype if self.dtype is None else self.dtype.__name__))

        if not self.index.is_unique:
            raise ValueError(
                "This MSA's index labels must be unique.")
        if not other.index.is_unique:
            raise ValueError(
                "`other`'s index labels must be unique.")

    def _get_join_index(self, other, how):
        if how == 'strict':
            diff = self.index.sym_diff(other.index)
            if len(diff) > 0:
                raise ValueError(
                    "Index labels must all match with `how='strict'`")

            diff = self.positional_metadata.columns.sym_diff(
                other.positional_metadata.columns)

            if not self.has_positional_metadata():
                del self.positional_metadata
            if not other.has_positional_metadata():
                del other.positional_metadata

            if len(diff) > 0:
                raise ValueError(
                    "Positional metadata columns must all match with "
                    "`how='strict'`")

            join_index = self.index
            concat_kwargs = {'join': 'inner'}
        elif how == 'inner':
            join_index = self.index.intersection(other.index)
            concat_kwargs = {'join': 'inner'}
        elif how == 'outer':
            join_index = self.index.union(other.index)
            concat_kwargs = {'join': 'outer'}
        elif how == 'left':
            join_index = self.index
            concat_kwargs = {'join_axes': [self.positional_metadata.columns]}
        else:  # how='right'
            join_index = other.index
            concat_kwargs = {'join_axes': [other.positional_metadata.columns]}

        return join_index, concat_kwargs

    def _get_sequence_for_join(self, label):
        if label in self.index:
            # TODO: use .loc when it is implemented.
            return self._get_sequence(self.index.get_loc(label))
        else:
            return self.dtype(
                self.dtype.default_gap_char * self.shape.position)

    def sort(self, level=None, ascending=True):
        """Sort sequences by index label in-place.

        Parameters
        ----------
        level : int or object, optional
            Index level to sort on when index is a ``pd.MultiIndex``. Does
            nothing otherwise.
        ascending: bool, optional
            If ``False``, sort in descending (i.e., reverse) order.

        See Also
        --------
        index
        reassign_index
        pandas.Series.sort_index

        Notes
        -----
        This is a passthrough to ``pd.Series.sort_index`` internally.

        Examples
        --------
        Create a ``TabularMSA`` object:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'c'}),
        ...         DNA('AC-', metadata={'id': 'b'}),
        ...         DNA('AC-', metadata={'id': 'a'})]
        >>> msa = TabularMSA(seqs, minter='id')

        Sort the sequences in alphabetical order by sequence identifier:

        >>> msa.sort()
        >>> msa == TabularMSA([DNA('AC-', metadata={'id': 'a'}),
        ...                    DNA('AC-', metadata={'id': 'b'}),
        ...                    DNA('ACG', metadata={'id': 'c'})], minter='id')
        True

        Note that since the sort is in-place, the ``TabularMSA`` object is
        modified (a new object is **not** returned).

        """
        series = self._seqs.sort_index(ascending=ascending, level=level)
        self._seqs = series

    @experimental(as_of='0.4.0-dev')
    def to_dict(self):
        """Create a ``dict`` from this ``TabularMSA``.

        Returns
        -------
        dict
            Dictionary constructed from the index labels and sequences in this
            ``TabularMSA``.

        Raises
        ------
        ValueError
            If index labels are not unique.

        See Also
        --------
        from_dict
        index
        reassign_index

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACGT'), DNA('A--T')]
        >>> msa = TabularMSA(seqs, index=['a', 'b'])
        >>> dictionary = msa.to_dict()
        >>> dictionary == {'a': DNA('ACGT'), 'b': DNA('A--T')}
        True

        """
        if self.index.is_unique:
            return self._seqs.to_dict()
        else:
            raise ValueError("Cannot convert to dict. Index labels are not"
                             " unique.")

    def _get_sequence(self, i):
        return self._seqs.iloc[i]

    def _get_position(self, i):
        seq = Sequence.concat([s[i] for s in self._seqs], how='outer')
        if self.has_positional_metadata():
            seq.metadata = dict(self.positional_metadata.iloc[i])
        return seq

    def _is_sequence_axis(self, axis):
        if axis == 'sequence' or axis == 0:
            return True
        elif axis == 'position' or axis == 1:
            return False
        else:
            raise ValueError(
                "`axis` must be 'sequence' (0) or 'position' (1), not %r"
                % axis)

    @overrides(PositionalMetadataMixin)
    def _positional_metadata_axis_len_(self):
        return self.shape.position
