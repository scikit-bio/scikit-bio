# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import copy

import numpy as np
import pandas as pd
import scipy.stats

from skbio._base import SkbioObject
from skbio.metadata._mixin import MetadataMixin, PositionalMetadataMixin
from skbio.sequence import Sequence
from skbio.sequence._grammared_sequence import GrammaredSequence
from skbio.util._decorator import experimental, classonlymethod, overrides
from skbio.util._misc import resolve_key
from skbio.alignment._indexing import TabularMSAILoc, TabularMSALoc

from skbio.alignment._repr import _TabularMSAReprBuilder


_Shape = collections.namedtuple('Shape', ['sequence', 'position'])


class TabularMSA(MetadataMixin, PositionalMetadataMixin, SkbioObject):
    """Store a multiple sequence alignment in tabular (row/column) form.

    Parameters
    ----------
    sequences : iterable of GrammaredSequence, TabularMSA
        Aligned sequences in the MSA. Sequences must all be the same type and
        length. For example, `sequences` could be an iterable of ``DNA``,
        ``RNA``, or ``Protein`` sequences. If `sequences` is a ``TabularMSA``,
        its `metadata`, `positional_metadata`, and `index` will be used unless
        overridden by parameters `metadata`, `positional_metadata`, and
        `minter`/`index`, respectively.
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
        a key into each sequence's ``metadata`` attribute. Note that `minter`
        cannot be combined with `index`.
    index : pd.Index consumable, optional
        Index containing labels for `sequences`. Must be the same length as
        `sequences`. Must be able to be passed directly to ``pd.Index``
        constructor. Note that `index` cannot be combined with `minter` and the
        contents of `index` must be hashable.

    Raises
    ------
    ValueError
        If `minter` and `index` are both provided.
    ValueError
        If `index` is not the same length as `sequences`.
    TypeError
        If `sequences` contains an object that isn't a ``GrammaredSequence``.
    TypeError
        If `sequences` does not contain exactly the same type of
        ``GrammaredSequence`` objects.
    ValueError
        If `sequences` does not contain ``GrammaredSequence`` objects of the
        same length.

    See Also
    --------
    skbio.sequence.DNA
    skbio.sequence.RNA
    skbio.sequence.Protein
    pandas.DataFrame
    pandas.Index
    reassign_index

    Notes
    -----
    If neither `minter` nor `index` are provided, default index labels will be
    used: ``pd.RangeIndex(start=0, stop=len(sequences), step=1)``.

    Examples
    --------
    Create a ``TabularMSA`` object with three DNA sequences and four positions:

    >>> from skbio import DNA, TabularMSA
    >>> seqs = [
    ...     DNA('ACGT'),
    ...     DNA('AG-T'),
    ...     DNA('-C-T')
    ... ]
    >>> msa = TabularMSA(seqs)
    >>> msa
    TabularMSA[DNA]
    ---------------------
    Stats:
        sequence count: 3
        position count: 4
    ---------------------
    ACGT
    AG-T
    -C-T

    Since `minter` or `index` wasn't provided, the MSA has default index
    labels:

    >>> msa.index
    RangeIndex(start=0, stop=3, step=1)

    Create an MSA with metadata, positional metadata, and non-default index
    labels:

    >>> msa = TabularMSA(seqs, index=['seq1', 'seq2', 'seq3'],
    ...                  metadata={'id': 'msa-id'},
    ...                  positional_metadata={'prob': [3, 4, 2, 2]})
    >>> msa
    TabularMSA[DNA]
    --------------------------
    Metadata:
        'id': 'msa-id'
    Positional metadata:
        'prob': <dtype: int64>
    Stats:
        sequence count: 3
        position count: 4
    --------------------------
    ACGT
    AG-T
    -C-T
    >>> msa.index
    Index(['seq1', 'seq2', 'seq3'], dtype='object')

    """
    default_write_format = 'fasta'
    __hash__ = None

    @property
    @experimental(as_of='0.4.1')
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
        return type(self._get_sequence_iloc_(0)) if len(self) > 0 else None

    @property
    @experimental(as_of='0.4.1')
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
            position_count = len(self._get_sequence_iloc_(0))
        else:
            position_count = 0

        return _Shape(sequence=sequence_count, position=position_count)

    @property
    @experimental(as_of='0.4.1')
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
        ...         DNA('AC-', metadata={'id': 'b'}),
        ...         DNA('AC-', metadata={'id': 'c'})]
        >>> msa = TabularMSA(seqs, minter='id')

        Retrieve index:

        >>> msa.index
        Index(['a', 'b', 'c'], dtype='object')

        Set index:

        >>> msa.index = ['seq1', 'seq2', 'seq3']
        >>> msa.index
        Index(['seq1', 'seq2', 'seq3'], dtype='object')

        Deleting the index resets it to the ``TabularMSA`` constructor's
        default:

        >>> del msa.index
        >>> msa.index
        RangeIndex(start=0, stop=3, step=1)

        """
        return self._seqs.index

    @index.setter
    def index(self, index):
        # Cast to Index to identify tuples as a MultiIndex to match
        # pandas constructor. Just setting would make an index of tuples.
        if not isinstance(index, pd.Index):
            index = pd.Index(index)
        self._seqs.index = index

    @index.deleter
    def index(self):
        # Create a memory-efficient integer index as the default MSA index.
        self._seqs.index = pd.RangeIndex(start=0, stop=len(self), step=1)

    @property
    @experimental(as_of="0.4.1")
    def loc(self):
        """Slice the MSA on first axis by index label, second axis by position.

        This will return an object with the following interface:

        .. code-block:: python

           msa.loc[seq_idx]
           msa.loc[seq_idx, pos_idx]
           msa.loc(axis='sequence')[seq_idx]
           msa.loc(axis='position')[pos_idx]

        Parameters
        ----------
        seq_idx : label, slice, 1D array_like (bool or label)
            Slice the first axis of the MSA. When this value is a scalar, a
            sequence of ``msa.dtype`` will be returned. This may be further
            sliced by `pos_idx`.
        pos_idx : (same as seq_idx), optional
            Slice the second axis of the MSA. When this value is a scalar, a
            sequence of type :class:`skbio.sequence.Sequence` will be returned.
            This represents a column of the MSA and may have been additionally
            sliced by `seq_idx`.
        axis : {'sequence', 'position', 0, 1, None}, optional
            Limit the axis to slice on. When set, a tuple as the argument will
            no longer be split into `seq_idx` and `pos_idx`.

        Returns
        -------
        TabularMSA, GrammaredSequence, Sequence
            A ``TabularMSA`` is returned when `seq_idx` and `pos_idx` are
            non-scalars. A ``GrammaredSequence`` of type ``msa.dtype`` is
            returned when `seq_idx` is a scalar (this object will match the
            dtype of the MSA). A ``Sequence`` is returned when `seq_idx` is
            non-scalar and `pos_idx` is scalar.

        See Also
        --------
        iloc
        __getitem__

        Notes
        -----
        If the slice operation results in a ``TabularMSA`` without any
        sequences, the MSA's ``positional_metadata`` will be unset.

        When the MSA's index is a ``pd.MultiIndex`` a tuple may be given to
        `seq_idx` to indicate the slicing operations to perform on each
        component index.

        Examples
        --------
        First we need to set up an MSA to slice:

        >>> from skbio import TabularMSA, DNA
        >>> msa = TabularMSA([DNA("ACGT"), DNA("A-GT"), DNA("AC-T"),
        ...                   DNA("ACGA")], index=['a', 'b', 'c', 'd'])
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 4
            position count: 4
        ---------------------
        ACGT
        A-GT
        AC-T
        ACGA
        >>> msa.index
        Index(['a', 'b', 'c', 'd'], dtype='object')


        When we slice by a scalar we get the original sequence back out of the
        MSA:

        >>> msa.loc['b']
        DNA
        --------------------------
        Stats:
            length: 4
            has gaps: True
            has degenerates: False
            has definites: True
            GC-content: 33.33%
        --------------------------
        0 A-GT

        Similarly when we slice the second axis by a scalar we get a column of
        the MSA:

        >>> msa.loc[..., 1]
        Sequence
        -------------
        Stats:
            length: 4
        -------------
        0 C-CC

        Note: we return an ``skbio.Sequence`` object because the column of an
        alignment has no biological meaning and many operations defined for the
        MSA's sequence `dtype` would be meaningless.

        When we slice both axes by a scalar, operations are applied left to
        right:

        >>> msa.loc['a', 0]
        DNA
        --------------------------
        Stats:
            length: 1
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 0.00%
        --------------------------
        0 A

        In other words, it exactly matches slicing the resulting sequence
        object directly:

        >>> msa.loc['a'][0]
        DNA
        --------------------------
        Stats:
            length: 1
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 0.00%
        --------------------------
        0 A

        When our slice is non-scalar we get back an MSA of the same `dtype`:

        >>> msa.loc[['a', 'c']]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        ACGT
        AC-T

        We can similarly slice out a column of that:

        >>> msa.loc[['a', 'c'], 2]
        Sequence
        -------------
        Stats:
            length: 2
        -------------
        0 G-

        Slice syntax works as well:

        >>> msa.loc[:'c']
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 4
        ---------------------
        ACGT
        A-GT
        AC-T

        Notice how the end label is included in the results. This is different
        from how positional slices behave:

        >>> msa.loc[[True, False, False, True], 2:3]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 1
        ---------------------
        G
        G

        Here we sliced the first axis by a boolean vector, but then restricted
        the columns to a single column. Because the second axis was given a
        nonscalar we still recieve an MSA even though only one column is
        present.

        Duplicate labels can be an unfortunate reality in the real world,
        however `loc` is capable of handling this:

        >>> msa.index = ['a', 'a', 'b', 'c']

        Notice how the label 'a' happens twice. If we were to access 'a' we get
        back an MSA with both sequences:

        >>> msa.loc['a']
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        ACGT
        A-GT

        Remember that `iloc` can always be used to differentiate sequences with
        duplicate labels.

        More advanced slicing patterns are possible with different index types.

        Let's use a `pd.MultiIndex`:

        >>> msa.index = [('a', 0), ('a', 1), ('b', 0), ('b', 1)]

        Here we will explicitly set the axis that we are slicing by to make
        things easier to read:

        >>> msa.loc(axis='sequence')['a', 0]
        DNA
        --------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 50.00%
        --------------------------
        0 ACGT

        This selected the first sequence because the complete label was
        provided. In other words `('a', 0)` was treated as a scalar for this
        index.

        We can also slice along the component indices of the multi-index:

        >>> msa.loc(axis='sequence')[:, 1]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        A-GT
        ACGA

        If we were to do that again without the `axis` argument, it would look
        like this:

        >>> msa.loc[(slice(None), 1), ...]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        A-GT
        ACGA

        Notice how we needed to specify the second axis. If we had left that
        out we would have simply gotten the 2nd column back instead. We also
        lost the syntactic sugar for slice objects. These are a few of the
        reasons specifying the `axis` preemptively can be useful.

        """
        return self._loc

    @property
    @experimental(as_of="0.4.1")
    def iloc(self):
        """Slice the MSA on either axis by index position.

        This will return an object with the following interface:

        .. code-block:: python

           msa.iloc[seq_idx]
           msa.iloc[seq_idx, pos_idx]
           msa.iloc(axis='sequence')[seq_idx]
           msa.iloc(axis='position')[pos_idx]

        Parameters
        ----------
        seq_idx : int, slice, iterable (int and slice), 1D array_like (bool)
            Slice the first axis of the MSA. When this value is a scalar, a
            sequence of ``msa.dtype`` will be returned. This may be further
            sliced by `pos_idx`.
        pos_idx : (same as seq_idx), optional
            Slice the second axis of the MSA. When this value is a scalar, a
            sequence of type :class:`skbio.sequence.Sequence` will be returned.
            This represents a column of the MSA and may have been additionally
            sliced by `seq_idx`.
        axis : {'sequence', 'position', 0, 1, None}, optional
            Limit the axis to slice on. When set, a tuple as the argument will
            no longer be split into `seq_idx` and `pos_idx`.

        Returns
        -------
        TabularMSA, GrammaredSequence, Sequence
            A ``TabularMSA`` is returned when `seq_idx` and `pos_idx` are
            non-scalars. A ``GrammaredSequence`` of type ``msa.dtype`` is
            returned when `seq_idx` is a scalar (this object will match the
            dtype of the MSA). A ``Sequence`` is returned when `seq_idx` is
            non-scalar and `pos_idx` is scalar.

        See Also
        --------
        __getitem__
        loc

        Notes
        -----
        If the slice operation results in a ``TabularMSA`` without any
        sequences, the MSA's ``positional_metadata`` will be unset.

        Examples
        --------
        First we need to set up an MSA to slice:

        >>> from skbio import TabularMSA, DNA
        >>> msa = TabularMSA([DNA("ACGT"), DNA("A-GT"), DNA("AC-T"),
        ...                   DNA("ACGA")])
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 4
            position count: 4
        ---------------------
        ACGT
        A-GT
        AC-T
        ACGA

        When we slice by a scalar we get the original sequence back out of the
        MSA:

        >>> msa.iloc[1]
        DNA
        --------------------------
        Stats:
            length: 4
            has gaps: True
            has degenerates: False
            has definites: True
            GC-content: 33.33%
        --------------------------
        0 A-GT

        Similarly when we slice the second axis by a scalar we get a column of
        the MSA:

        >>> msa.iloc[..., 1]
        Sequence
        -------------
        Stats:
            length: 4
        -------------
        0 C-CC

        Note: we return an ``skbio.Sequence`` object because the column of an
        alignment has no biological meaning and many operations defined for the
        MSA's sequence `dtype` would be meaningless.

        When we slice both axes by a scalar, operations are applied left to
        right:

        >>> msa.iloc[0, 0]
        DNA
        --------------------------
        Stats:
            length: 1
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 0.00%
        --------------------------
        0 A

        In other words, it exactly matches slicing the resulting sequence
        object directly:

        >>> msa.iloc[0][0]
        DNA
        --------------------------
        Stats:
            length: 1
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 0.00%
        --------------------------
        0 A

        When our slice is non-scalar we get back an MSA of the same `dtype`:

        >>> msa.iloc[[0, 2]]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        ACGT
        AC-T

        We can similarly slice out a column of that:

        >>> msa.iloc[[0, 2], 2]
        Sequence
        -------------
        Stats:
            length: 2
        -------------
        0 G-

        Slice syntax works as well:

        >>> msa.iloc[:3]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 4
        ---------------------
        ACGT
        A-GT
        AC-T

        We can also use boolean vectors:

        >>> msa.iloc[[True, False, False, True], 2:3]
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 1
        ---------------------
        G
        G

        Here we sliced the first axis by a boolean vector, but then restricted
        the columns to a single column. Because the second axis was given a
        nonscalar we still recieve an MSA even though only one column is
        present.

        """
        return self._iloc

    @classonlymethod
    @experimental(as_of="0.4.1")
    def from_dict(cls, dictionary):
        """Create a ``TabularMSA`` from a ``dict``.

        Parameters
        ----------
        dictionary : dict
            Dictionary mapping keys to ``GrammaredSequence`` sequence objects.
            The ``TabularMSA`` object will have its index labels set
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
        >>> msa.shape
        Shape(sequence=2, position=4)
        >>> 'a' in msa
        True
        >>> 'b' in msa
        True

        """
        # Python 3 guarantees same order of iteration as long as no
        # modifications are made to the dictionary between calls:
        #     https://docs.python.org/3/library/stdtypes.html#
        #         dictionary-view-objects
        return cls(dictionary.values(), index=dictionary.keys())

    @experimental(as_of='0.4.1')
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

        # Give a better error message than the one raised by `extend` (it
        # references `reset_index`, which isn't a constructor parameter).
        if minter is not None and index is not None:
            raise ValueError(
                "Cannot use both `minter` and `index` at the same time.")
        self._seqs = pd.Series([])
        self.extend(sequences, minter=minter, index=index,
                    reset_index=minter is None and index is None)

        MetadataMixin._init_(self, metadata=metadata)
        PositionalMetadataMixin._init_(
            self, positional_metadata=positional_metadata)

        # Set up our indexers
        self._loc = TabularMSALoc(self)
        self._iloc = TabularMSAILoc(self)

    def _constructor_(self, sequences=NotImplemented, metadata=NotImplemented,
                      positional_metadata=NotImplemented,
                      index=NotImplemented):
        """Return new copy of the MSA with overridden properties.

        NotImplemented is used as a sentinel so that None may be used to
        override values.
        """
        if metadata is NotImplemented:
            if self.has_metadata():
                metadata = self.metadata
            else:
                metadata = None

        if positional_metadata is NotImplemented:
            if self.has_positional_metadata():
                positional_metadata = self.positional_metadata
            else:
                positional_metadata = None

        if index is NotImplemented:
            if isinstance(sequences, pd.Series):
                index = sequences.index
            else:
                index = self.index

        if sequences is NotImplemented:
            sequences = self._seqs

        sequences = [copy.copy(s) for s in sequences]

        return self.__class__(sequences, metadata=metadata,
                              positional_metadata=positional_metadata,
                              index=index)

    @experimental(as_of='0.4.1')
    def __repr__(self):
        """String summary of this MSA."""
        pep8_line_length_limit = 79
        length_taken_by_docstring_indent = 8
        width = pep8_line_length_limit - length_taken_by_docstring_indent
        return _TabularMSAReprBuilder(
            msa=self,
            width=width,
            indent=4).build()

    def _repr_stats(self):
        return [("sequence count", str(self.shape.sequence)),
                ("position count", str(self.shape.position))]

    @experimental(as_of='0.4.1')
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
        # TODO: change for #1198
        return self.shape.position > 0

    @experimental(as_of='0.4.1')
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

    @experimental(as_of='0.4.1')
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

    @experimental(as_of='0.4.1')
    def __iter__(self):
        """Iterate over sequences in the MSA.

        Yields
        ------
        GrammaredSequence
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

    @experimental(as_of='0.4.1')
    def __reversed__(self):
        """Iterate in reverse order over sequences in the MSA.

        Yields
        ------
        GrammaredSequence
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

    @experimental(as_of='0.4.1')
    def __str__(self):
        """String summary of this MSA."""
        return self.__repr__()

    @experimental(as_of='0.4.1')
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

    @experimental(as_of='0.4.1')
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

    @experimental(as_of='0.4.1')
    def __copy__(self):
        """Return a shallow copy of this MSA.

        Returns
        -------
        TabularMSA
            Shallow copy of this MSA. Sequence objects will be shallow-copied.

        See Also
        --------
        __deepcopy__

        Examples
        --------
        >>> import copy
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa_copy = copy.copy(msa)
        >>> msa_copy == msa
        True
        >>> msa_copy is msa
        False

        """
        msa_copy = self._constructor_()

        msa_copy._metadata = MetadataMixin._copy_(self)
        msa_copy._positional_metadata = PositionalMetadataMixin._copy_(self)

        return msa_copy

    @experimental(as_of='0.4.1')
    def __deepcopy__(self, memo):
        """Return a deep copy of this MSA.

        Returns
        -------
        TabularMSA
            Deep copy of this MSA. Sequence objects will be deep-copied.

        See Also
        --------
        __copy__

        Examples
        --------
        >>> import copy
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa_copy = copy.deepcopy(msa)
        >>> msa_copy == msa
        True
        >>> msa_copy is msa
        False

        """
        seqs = (copy.deepcopy(seq, memo) for seq in self._seqs)
        msa_copy = self._constructor_(sequences=seqs)

        msa_copy._metadata = MetadataMixin._deepcopy_(self, memo)
        msa_copy._positional_metadata = \
            PositionalMetadataMixin._deepcopy_(self, memo)

        return msa_copy

    @experimental(as_of="0.4.1")
    def __getitem__(self, indexable):
        """Slice the MSA on either axis.

        This is a pass-through for :func:`skbio.alignment.TabularMSA.iloc`.
        Please refer to the associated documentation.

        See Also
        --------
        iloc
        loc

        Notes
        -----
        Axis restriction is not possible for this method.

        To slice by labels, use ``loc``.

        """
        return self.iloc[indexable]

    # Helpers for TabularMSAILoc and TabularMSALoc
    def _get_sequence_iloc_(self, i):
        return self._seqs.iloc[i]

    def _slice_sequences_iloc_(self, i):
        new_seqs = self._seqs.iloc[i]
        # TODO: change for #1198
        if len(new_seqs) == 0:
            return self._constructor_(new_seqs, positional_metadata=None)
        return self._constructor_(new_seqs)

    def _get_sequence_loc_(self, ids):
        new_seqs = self._seqs.loc[ids]
        if type(new_seqs) is self.dtype:
            return new_seqs
        else:
            # Thanks CategoricalIndex, you understand no such thing as a scalar
            if len(new_seqs) == 1:
                return new_seqs.iloc[0]
            else:
                # This was a common failure mode; shouldn't happen anymore, but
                # it could strike again.
                raise AssertionError(
                    "Something went wrong with the index %r provided to"
                    " `_get_sequence_loc_`, please report this stack trace to"
                    "\nhttps://github.com/biocore/scikit-bio/issues" % ids)

    def _slice_sequences_loc_(self, ids):
        new_seqs = self._seqs.loc[ids]
        try:
            # TODO: change for #1198
            if len(new_seqs) == 0:
                return self._constructor_(new_seqs, positional_metadata=None)
            return self._constructor_(new_seqs)
        except TypeError:  # NaN hit the constructor, key was bad... probably
            raise KeyError("Part of `%r` was not in the index." % ids)

    def _get_position_(self, i, ignore_metadata=False):
        if ignore_metadata:
            return Sequence(''.join([str(s[i]) for s in self._seqs]))

        seq = Sequence.concat([s[i] for s in self._seqs], how='outer')
        # TODO: change for #1198
        if len(self) and self.has_positional_metadata():
            seq.metadata = dict(self.positional_metadata.iloc[i])
        return seq

    def _slice_positions_(self, i):
        seqs = self._seqs.apply(lambda seq: seq[i])
        # TODO: change for #1198
        pm = None
        if len(self) and self.has_positional_metadata():
            pm = self.positional_metadata.iloc[i]
        return self._constructor_(seqs, positional_metadata=pm)
    # end of helpers

    @experimental(as_of='0.4.1')
    def iter_positions(self, reverse=False, ignore_metadata=False):
        """Iterate over positions (columns) in the MSA.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, iterate over positions in reverse order.
        ignore_metadata : bool, optional
            If ``True``, ``Sequence.metadata`` and
            ``Sequence.positional_metadata`` will not be included. This can
            significantly improve performance if metadata is not needed.

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
        metadata stored as ``metadata`` unless ``ignore_metadata`` is set to
        ``True``.

        Sequences will have their positional metadata concatenated using an
        outer join unless ``ignore_metadata`` is set to ``True``. See
        ``Sequence.concat(how='outer')`` for details.

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

        return (self._get_position_(index, ignore_metadata=ignore_metadata)
                for index in indices)

    @experimental(as_of='0.4.1')
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
        --------------------------
        Positional metadata:
            'prob': <dtype: int64>
        Stats:
            length: 5
            has gaps: True
            has degenerates: False
            has definites: True
            GC-content: 33.33%
        --------------------------
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
        for position in self.iter_positions(ignore_metadata=True):
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

    def _build_inverse_shannon_uncertainty_f(self, include_gaps):
        base = len(self.dtype.definite_chars)
        if include_gaps:
            # Increment the base by one to reflect the possible inclusion of
            # the default gap character.
            base += 1

        def f(p):
            freqs = list(p.frequencies().values())
            return 1. - scipy.stats.entropy(freqs, base=base)
        return f

    @experimental(as_of='0.4.1')
    def conservation(self, metric='inverse_shannon_uncertainty',
                     degenerate_mode='error', gap_mode='nan'):
        """Apply metric to compute conservation for all alignment positions

        Parameters
        ----------
        metric : {'inverse_shannon_uncertainty'}, optional
            Metric that should be applied for computing conservation. Resulting
            values should be larger when a position is more conserved.
        degenerate_mode : {'nan', 'error'}, optional
            Mode for handling positions with degenerate characters. If
            ``"nan"``, positions with degenerate characters will be assigned a
            conservation score of ``np.nan``. If ``"error"``, an
            error will be raised if one or more degenerate characters are
            present.
        gap_mode : {'nan', 'ignore', 'error', 'include'}, optional
            Mode for handling positions with gap characters. If ``"nan"``,
            positions with gaps will be assigned a conservation score of
            ``np.nan``. If ``"ignore"``, positions with gaps will be filtered
            to remove gaps before ``metric`` is applied. If ``"error"``, an
            error will be raised if one or more gap characters are present. If
            ``"include"``, conservation will be computed on alignment positions
            with gaps included. In this case, it is up to the metric to ensure
            that gaps are handled as they should be or to raise an error if
            gaps are not supported by that metric.

        Returns
        -------
        np.array of floats
            Values resulting from the application of ``metric`` to each
            position in the alignment.

        Raises
        ------
        ValueError
            If an unknown ``metric``, ``degenerate_mode`` or ``gap_mode`` is
            provided.
        ValueError
            If any degenerate characters are present in the alignment when
            ``degenerate_mode`` is ``"error"``.
        ValueError
            If any gaps are present in the alignment when ``gap_mode`` is
            ``"error"``.

        Notes
        -----
        Users should be careful interpreting results when
        ``gap_mode = "include"`` as the results may be misleading. For example,
        as pointed out in [1]_, a protein alignment position composed of 90%
        gaps and 10% tryptophans would score as more highly conserved than a
        position composed of alanine and glycine in equal frequencies with the
        ``"inverse_shannon_uncertainty"`` metric.

        ``gap_mode = "include"`` will result in all gap characters being
        recoded to ``TabularMSA.dtype.default_gap_char``. Because no
        conservation metrics that we are aware of consider different gap
        characters differently (e.g., none of the metrics described in [1]_),
        they are all treated the same within this method.

        The ``inverse_shannon_uncertainty`` metric is simply one minus
        Shannon's uncertainty metric. This method uses the inverse of Shannon's
        uncertainty so that larger values imply higher conservation. Shannon's
        uncertainty is also referred to as Shannon's entropy, but when making
        computations from symbols, as is done here, "uncertainty" is the
        preferred term ([2]_).

        References
        ----------
        .. [1] Valdar WS. Scoring residue conservation. Proteins. (2002)
        .. [2] Schneider T. Pitfalls in information theory (website, ca. 2015).
           https://schneider.ncifcrf.gov/glossary.html#Shannon_entropy

        """

        if gap_mode not in {'nan', 'error', 'include', 'ignore'}:
            raise ValueError("Unknown gap_mode provided: %s" % gap_mode)

        if degenerate_mode not in {'nan', 'error'}:
            raise ValueError("Unknown degenerate_mode provided: %s" %
                             degenerate_mode)

        if metric not in {'inverse_shannon_uncertainty'}:
            raise ValueError("Unknown metric provided: %s" %
                             metric)

        if self.shape[0] == 0:
            # handle empty alignment to avoid error on lookup of character sets
            return np.array([])

        # Since the only currently allowed metric is
        # inverse_shannon_uncertainty, and we already know that a valid metric
        # was provided, we just define metric_f here. When additional metrics
        # are supported, this will be handled differently (e.g., via a lookup
        # or if/elif/else).
        metric_f = self._build_inverse_shannon_uncertainty_f(
                        gap_mode == 'include')

        result = []
        for p in self.iter_positions(ignore_metadata=True):
            cons = None
            # cast p to self.dtype for access to gap/degenerate related
            # functionality
            pos_seq = self.dtype(p)

            # handle degenerate characters if present
            if pos_seq.has_degenerates():
                if degenerate_mode == 'nan':
                    cons = np.nan
                else:  # degenerate_mode == 'error' is the only choice left
                    degenerate_chars = pos_seq[pos_seq.degenerates()]
                    raise ValueError("Conservation is undefined for positions "
                                     "with degenerate characters. The "
                                     "following degenerate characters were "
                                     "observed: %s." % degenerate_chars)

            # handle gap characters if present
            if pos_seq.has_gaps():
                if gap_mode == 'nan':
                    cons = np.nan
                elif gap_mode == 'error':
                    raise ValueError("Gap characters present in alignment.")
                elif gap_mode == 'ignore':
                    pos_seq = pos_seq.degap()
                else:  # gap_mode == 'include' is the only choice left
                    # Recode all gap characters with pos_seq.default_gap_char.
                    pos_seq = pos_seq.replace(pos_seq.gaps(),
                                              pos_seq.default_gap_char)

            if cons is None:
                cons = metric_f(pos_seq)

            result.append(cons)

        return np.array(result)

    @experimental(as_of='0.4.1')
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
            seq_iterator = self.iter_positions(ignore_metadata=True)
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
            gap_freqs.append(sum(freqs.values()))

        gap_freqs = np.asarray(gap_freqs, dtype=float if relative else int)

        if relative:
            gap_freqs /= length

        return gap_freqs

    @experimental(as_of='0.4.1')
    def reassign_index(self, mapping=None, minter=None):
        """Reassign index labels to sequences in this MSA.

        Parameters
        ----------
        mapping : dict or callable, optional
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
        If neither `mapping` nor `minter` are provided, index labels will be
        reset to the ``TabularMSA`` constructor's default.

        Examples
        --------
        Create a ``TabularMSA`` object with default index labels:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'a'}),
        ...         DNA('AC-', metadata={'id': 'b'}),
        ...         DNA('CCG', metadata={'id': 'c'})]
        >>> msa = TabularMSA(seqs)
        >>> msa.index
        RangeIndex(start=0, stop=3, step=1)

        Assign new index to the MSA using each sequence's ID as a label:

        >>> msa.reassign_index(minter='id')
        >>> msa.index
        Index(['a', 'b', 'c'], dtype='object')

        Assign default index:

        >>> msa.reassign_index()
        >>> msa.index
        RangeIndex(start=0, stop=3, step=1)

        Alternatively, a mapping of existing labels to new labels may be passed
        via `mapping`:

        >>> msa.reassign_index(mapping={0: 'seq1', 1: 'seq2'})
        >>> msa.index
        Index(['seq1', 'seq2', 2], dtype='object')

        """
        if mapping is not None and minter is not None:
            raise ValueError(
                "Cannot use both `mapping` and `minter` at the same time.")

        if mapping is not None:
            if isinstance(mapping, dict):
                self.index = [mapping[label] if label in mapping else label
                              for label in self.index]
            elif callable(mapping):
                self.index = [mapping(label) for label in self.index]
            else:
                raise TypeError(
                    "`mapping` must be a dict or callable, not type %r"
                    % type(mapping).__name__)
        elif minter is not None:
            self.index = [resolve_key(seq, minter) for seq in self._seqs]
        else:
            del self.index

    @experimental(as_of='0.4.1')
    def append(self, sequence, minter=None, index=None, reset_index=False):
        """Append a sequence to the MSA without recomputing alignment.

        Parameters
        ----------
        sequence : GrammaredSequence
            Sequence to be appended. Must match the dtype of the MSA and the
            number of positions in the MSA.
        minter : callable or metadata key, optional
            Used to create an index label for the sequence being appended. If
            callable, it generates a label directly. Otherwise it's treated as
            a key into the sequence metadata. Note that `minter` cannot be
            combined with `index` nor `reset_index`.
        index : object, optional
            Index label to use for the appended sequence. Note that `index`
            cannot be combined with `minter` nor `reset_index`.
        reset_index : bool, optional
            If ``True``, this MSA's index is reset to the ``TabularMSA``
            constructor's default after appending. Note that `reset_index`
            cannot be combined with `minter` nor `index`.

        Raises
        ------
        ValueError
            If exactly one choice of `minter`, `index`, or `reset_index` is not
            provided.
        TypeError
            If the sequence object isn't a ``GrammaredSequence``.
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
        The MSA is not automatically re-aligned when a sequence is appended.
        Therefore, this operation is not necessarily meaningful on its own.

        Examples
        --------
        Create an MSA with a single sequence labeled ``'seq1'``:

        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACGT')], index=['seq1'])
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 1
            position count: 4
        ---------------------
        ACGT
        >>> msa.index
        Index(['seq1'], dtype='object')

        Append a new sequence to the MSA, providing its index label via
        `index`:

        >>> msa.append(DNA('AG-T'), index='seq2')
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 4
        ---------------------
        ACGT
        AG-T
        >>> msa.index
        Index(['seq1', 'seq2'], dtype='object')

        Append another sequence, this time resetting the MSA's index labels to
        the default with `reset_index`. Note that since the MSA's index is
        reset, we do not need to provide an index label for the new sequence
        via `index` or `minter`:

        >>> msa.append(DNA('ACGA'), reset_index=True)
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 4
        ---------------------
        ACGT
        AG-T
        ACGA
        >>> msa.index
        RangeIndex(start=0, stop=3, step=1)

        """
        if index is not None:
            index = [index]
        self.extend([sequence], minter=minter, index=index,
                    reset_index=reset_index)

    @experimental(as_of='0.4.1')
    def extend(self, sequences, minter=None, index=None, reset_index=False):
        """Extend this MSA with sequences without recomputing alignment.

        Parameters
        ----------
        sequences : iterable of GrammaredSequence
            Sequences to be appended. Must match the dtype of the MSA and the
            number of positions in the MSA.
        minter : callable or metadata key, optional
            Used to create index labels for the sequences being appended. If
            callable, it generates a label directly. Otherwise it's treated as
            a key into the sequence metadata. Note that `minter` cannot be
            combined with `index` nor `reset_index`.
        index : pd.Index consumable, optional
            Index labels to use for the appended sequences. Must be the same
            length as `sequences`. Must be able to be passed directly to
            ``pd.Index`` constructor. Note that `index` cannot be combined
            with `minter` nor `reset_index`.
        reset_index : bool, optional
            If ``True``, this MSA's index is reset to the ``TabularMSA``
            constructor's default after extending. Note that `reset_index`
            cannot be combined with `minter` nor `index`.

        Raises
        ------
        ValueError
            If exactly one choice of `minter`, `index`, or `reset_index` is not
            provided.
        ValueError
            If `index` is not the same length as `sequences`.
        TypeError
            If `sequences` contains an object that isn't a
            ``GrammaredSequence``.
        TypeError
            If `sequences` contains a type that does not match the dtype of the
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
        The MSA is not automatically re-aligned when appending sequences.
        Therefore, this operation is not necessarily meaningful on its own.

        Examples
        --------
        Create an MSA with a single sequence labeled ``'seq1'``:

        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACGT')], index=['seq1'])
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 1
            position count: 4
        ---------------------
        ACGT
        >>> msa.index
        Index(['seq1'], dtype='object')

        Extend the MSA with sequences, providing their index labels via
        `index`:

        >>> msa.extend([DNA('AG-T'), DNA('-G-T')], index=['seq2', 'seq3'])
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 4
        ---------------------
        ACGT
        AG-T
        -G-T
        >>> msa.index
        Index(['seq1', 'seq2', 'seq3'], dtype='object')

        Extend with more sequences, this time resetting the MSA's index labels
        to the default with `reset_index`. Note that since the MSA's index is
        reset, we do not need to provide index labels for the new sequences via
        `index` or `minter`:

        >>> msa.extend([DNA('ACGA'), DNA('AC-T'), DNA('----')],
        ...            reset_index=True)
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 6
            position count: 4
        ---------------------
        ACGT
        AG-T
        ...
        AC-T
        ----
        >>> msa.index
        RangeIndex(start=0, stop=6, step=1)

        """
        if sum([minter is not None,
                index is not None,
                bool(reset_index)]) != 1:
            raise ValueError(
                "Must provide exactly one of the following parameters: "
                "`minter`, `index`, `reset_index`")

        # Verify `sequences` first because `minter` could interact with each
        # sequence's `metadata`.
        sequences = list(sequences)
        self._assert_valid_sequences(sequences)

        if minter is not None:
            # Convert to Index to identify tuples as a MultiIndex instead of an
            # index of tuples.
            index = pd.Index([resolve_key(seq, minter) for seq in sequences])
        elif index is not None:
            # Convert to Index to identify tuples as a MultiIndex instead of an
            # index of tuples.
            if not isinstance(index, pd.Index):
                index = pd.Index(index)

            # pandas doesn't give a user-friendly error message if we pass
            # through.
            if len(sequences) != len(index):
                raise ValueError(
                    "Number of sequences (%d) must match index length (%d)" %
                    (len(sequences), len(index)))
        else:
            # Case for `reset_index=True`. We could simply set `index=None`
            # since it will be reset after appending below, but we can avoid a
            # memory spike if Series.append creates a new RangeIndex from
            # adjacent RangeIndexes in the future (pandas 0.18.0 creates an
            # Int64Index).
            index = pd.RangeIndex(start=len(self),
                                  stop=len(self) + len(sequences),
                                  step=1)

        if len(self):
            self._seqs = self._seqs.append(pd.Series(sequences, index=index))
        else:
            # Not using Series.append to avoid turning a RangeIndex supplied
            # via `index` parameter into an Int64Index (this happens in pandas
            # 0.18.0).
            self._seqs = pd.Series(sequences, index=index)

            # When extending a TabularMSA without sequences, the number of
            # positions in the TabularMSA may change from zero to non-zero. If
            # this happens, the TabularMSA's positional_metadata must be reset
            # to its default "empty" representation for the new number of
            # positions, otherwise the number of positions in the TabularMSA
            # and positional_metadata will differ.
            #
            # TODO: change for #1198
            if self.shape.position > 0:
                del self.positional_metadata

        if reset_index:
            self.reassign_index()

    def _assert_valid_sequences(self, sequences):
        if not sequences:
            return

        if len(self):
            expected_dtype = self.dtype
            expected_length = self.shape.position
        else:
            sequence = sequences[0]
            expected_dtype = type(sequence)
            if not issubclass(expected_dtype, GrammaredSequence):
                raise TypeError(
                    "Each sequence must be of type %r, not type %r"
                    % (GrammaredSequence.__name__, expected_dtype.__name__))
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
        .. note:: The following examples call `.sort()` on the joined MSA
           because there isn't a guaranteed ordering to the index. The joined
           MSA is sorted in these examples to make the output reproducible.
           When using this method with your own data, sorting the joined MSA is
           not necessary.

        Join MSAs by sequence:

        >>> from skbio import DNA, TabularMSA
        >>> msa1 = TabularMSA([DNA('AC'),
        ...                    DNA('A-')])
        >>> msa2 = TabularMSA([DNA('G-T'),
        ...                    DNA('T--')])
        >>> joined = msa1.join(msa2)
        >>> joined.sort()  # unnecessary in practice, see note above
        >>> joined
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 5
        ---------------------
        ACG-T
        A-T--

        Sequences are joined based on MSA index labels:

        >>> msa1 = TabularMSA([DNA('AC'),
        ...                    DNA('A-')], index=['a', 'b'])
        >>> msa2 = TabularMSA([DNA('G-T'),
        ...                    DNA('T--')], index=['b', 'a'])
        >>> joined = msa1.join(msa2)
        >>> joined.sort()  # unnecessary in practice, see note above
        >>> joined
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 2
            position count: 5
        ---------------------
        ACT--
        A-G-T
        >>> joined.index
        Index(['a', 'b'], dtype='object')

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
        >>> joined.sort()  # unnecessary in practice, see note above
        >>> joined
        TabularMSA[DNA]
        --------------------------
        Positional metadata:
            'col2': <dtype: int64>
        Stats:
            sequence count: 2
            position count: 5
        --------------------------
        ACT--
        A-G-T
        >>> joined.index
        Index(['a', 'b'], dtype='object')
        >>> joined.positional_metadata
           col2
        0     1
        1     2
        2     3
        3     4
        4     5

        When performing an outer join (``'outer'``, ``'left'``, or
        ``'right'``), unshared sequences are padded with gaps and unshared
        ``positional_metadata`` columns are padded with NaN:

        >>> joined = msa1.join(msa2, how='outer')
        >>> joined.sort()  # unnecessary in practice, see note above
        >>> joined
        TabularMSA[DNA]
        ----------------------------
        Positional metadata:
            'col1': <dtype: float64>
            'col2': <dtype: int64>
            'col3': <dtype: object>
        Stats:
            sequence count: 4
            position count: 5
        ----------------------------
        ACT--
        A-G-T
        -C---
        --ACG
        >>> joined.index
        Index(['a', 'b', 'c', 'z'], dtype='object')
        >>> joined.positional_metadata
           col1  col2 col3
        0  42.0     1  NaN
        1  43.0     2  NaN
        2   NaN     3    f
        3   NaN     4    o
        4   NaN     5    o

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
            if how == 'left':
                joined_positional_metadata = pd.concat(
                    [self.positional_metadata,
                     other.positional_metadata.reindex(
                        columns=self.positional_metadata.columns)],
                    ignore_index=True, sort=True)
            elif how == 'right':
                joined_positional_metadata = pd.concat(
                    [self.positional_metadata.reindex(
                        columns=other.positional_metadata.columns),
                     other.positional_metadata],
                    ignore_index=True, sort=True)
            else:
                joined_positional_metadata = pd.concat(
                    [self.positional_metadata, other.positional_metadata],
                    ignore_index=True, sort=True, **concat_kwargs)

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
            diff = self.index.symmetric_difference(other.index)
            if len(diff) > 0:
                raise ValueError(
                    "Index labels must all match with `how='strict'`")

            diff = self.positional_metadata.columns.symmetric_difference(
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
            return self.loc[label]
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
        Create a ``TabularMSA`` object with sequence identifiers as index
        labels:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'c'}),
        ...         DNA('AC-', metadata={'id': 'b'}),
        ...         DNA('AC-', metadata={'id': 'a'})]
        >>> msa = TabularMSA(seqs, minter='id')
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 3
        ---------------------
        ACG
        AC-
        AC-
        >>> msa.index
        Index(['c', 'b', 'a'], dtype='object')

        Sort the sequences in alphabetical order by index label:

        >>> msa.sort()
        >>> msa
        TabularMSA[DNA]
        ---------------------
        Stats:
            sequence count: 3
            position count: 3
        ---------------------
        AC-
        AC-
        ACG
        >>> msa.index
        Index(['a', 'b', 'c'], dtype='object')

        Note that since the sort is in-place, the ``TabularMSA`` object is
        modified (a new object is *not* returned).

        """
        series = self._seqs.sort_index(ascending=ascending, level=level)
        self._seqs = series
        positional_metadata = self.positional_metadata.sort_index(axis=1)
        self.positional_metadata = positional_metadata

    @experimental(as_of='0.4.1')
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
