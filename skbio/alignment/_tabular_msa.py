# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import collections

import numpy as np

from skbio._base import SkbioObject
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util import find_duplicates, OperationError, UniqueError
from skbio.util._decorator import experimental
from skbio.util._misc import resolve_key


_Shape = collections.namedtuple('Shape', ['sequence', 'position'])


class TabularMSA(SkbioObject):
    """Store a multiple sequence alignment in tabular (row/column) form.

    Parameters
    ----------
    sequences : iterable of alphabet-aware scikit-bio sequence objects
        Aligned sequences in the MSA. Sequences must be the same type, length,
        and have an alphabet. For example, `sequences` could be an iterable of
        ``DNA``, ``RNA``, or ``Protein`` objects.
    key : callable or metadata key, optional
        If provided, defines a unique, hashable key for each sequence in
        `sequences`. Can either be a callable accepting a single argument (each
        sequence) or a key into each sequence's ``metadata`` attribute.
    keys : iterable, optional
        An iterable of the same length as `sequences` containing unique,
        hashable elements. Each element will be used as the respective key for
        `sequences`.

    Raises
    ------
    ValueError
        If `key` and `keys` are both provided.
    ValueError
        If `keys` is not the same length as `sequences`.
    UniqueError
        If keys are not unique.

    See Also
    --------
    skbio.sequence.DNA
    skbio.sequence.RNA
    skbio.sequence.Protein

    Notes
    -----
    If `key` or `keys` are not provided, keys will not be set and certain
    operations requiring keys will raise an ``OperationError``.

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
        return self._dtype

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
        return self._shape

    @property
    @experimental(as_of='0.4.0-dev')
    def keys(self):
        """Keys in the order of sequences in the MSA.

        Raises
        ------
        OperationError
            If keys do not exist.

        See Also
        --------
        has_keys
        reindex

        Notes
        -----
        This property can be set and deleted. Keys are stored as an immutable
        ``numpy.ndarray``.

        Examples
        --------
        Create a ``TabularMSA`` object keyed by sequence identifier:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'a'}),
        ...         DNA('AC-', metadata={'id': 'b'})]
        >>> msa = TabularMSA(seqs, key='id')

        Retrieve keys:

        >>> msa.keys # doctest: +NORMALIZE_WHITESPACE
        array(['a', 'b'],
              dtype='<U1')

        Set keys:

        >>> msa.keys = ['seq1', 'seq2']
        >>> msa.keys # doctest: +NORMALIZE_WHITESPACE
        array(['seq1', 'seq2'],
              dtype='<U4')

        To make updates to a subset of the keys, first make a copy of the keys,
        update them, then set them again:

        >>> new_keys = msa.keys.copy()
        >>> new_keys[0] = 'new1'
        >>> msa.keys = new_keys
        >>> msa.keys # doctest: +NORMALIZE_WHITESPACE
        array(['new1', 'seq2'],
              dtype='<U4')

        Delete keys:

        >>> msa.has_keys()
        True
        >>> del msa.keys
        >>> msa.has_keys()
        False

        """
        if not self.has_keys():
            raise OperationError(
                "Keys do not exist. Use `reindex` to set them.")
        return self._keys

    @keys.setter
    def keys(self, keys):
        self.reindex(keys=keys)

    @keys.deleter
    def keys(self):
        self.reindex()

    @experimental(as_of='0.4.0-dev')
    def __init__(self, sequences, key=None, keys=None):
        sequences = iter(sequences)
        seq = next(sequences, None)

        dtype = None
        length = 0
        seqs = []
        if seq is not None:
            seqs.append(seq)
            dtype = type(seq)
            if not issubclass(dtype, IUPACSequence):
                raise TypeError(
                    "`sequences` must contain scikit-bio sequence objects "
                    "that have an alphabet, not type %r" % dtype.__name__)
            length = len(seq)

            for seq in sequences:
                if type(seq) is not dtype:
                    raise TypeError(
                        "`sequences` cannot contain mixed types. Type %r does "
                        "not match type %r" %
                        (type(seq).__name__, dtype.__name__))
                if len(seq) != length:
                    raise ValueError(
                        "`sequences` must contain sequences of the same "
                        "length: %r != %r" % (len(seq), length))
                seqs.append(seq)

        self._seqs = seqs
        self._dtype = dtype
        self._shape = _Shape(sequence=len(seqs), position=length)
        self.reindex(key=key, keys=keys)

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
        return self.shape.sequence

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
        # TODO implement me!
        return super(TabularMSA, self).__str__()

    @experimental(as_of='0.4.0-dev')
    def __eq__(self, other):
        """Determine if this MSA is equal to another.

        ``TabularMSA`` objects are equal if their sequences and keys are equal.

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
        >>> msa1 = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa1 == msa1
        True

        MSAs with different sequence characters are not equal:

        >>> msa2 = TabularMSA([DNA('ACG'), DNA('--G')])
        >>> msa1 == msa2
        False

        MSAs with different types of sequences (different ``dtype``) are not
        equal:

        >>> msa3 = TabularMSA([RNA('ACG'), RNA('AC-')])
        >>> msa1 == msa3
        False

        MSAs with different sequence metadata are not equal:

        >>> msa4 = TabularMSA([DNA('ACG', metadata={'id': 'a'}), DNA('AC-')])
        >>> msa1 == msa4
        False

        MSAs with different keys are not equal:

        >>> msa5 = TabularMSA([DNA('ACG'), DNA('AC-')], key=str)
        >>> msa1 == msa5
        False

        """
        if not isinstance(other, TabularMSA):
            return False

        # Use np.array_equal instead of (a == b).all():
        #   http://stackoverflow.com/a/10580782/3776794
        return ((self._seqs == other._seqs) and
                np.array_equal(self._keys, other._keys))

    @experimental(as_of='0.4.0-dev')
    def __ne__(self, other):
        """Determine if this MSA is not equal to another.

        ``TabularMSA`` objects are not equal if their sequences or keys are not
        equal.

        Parameters
        ----------
        other : TabularMSA
            MSA to test for inequality against.

        Returns
        -------
        bool
            Indicates whether this MSA is not equal to `other`.

        Examples
        --------
        >>> from skbio import DNA, RNA, TabularMSA
        >>> msa1 = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa1 != msa1
        False

        MSAs with different sequence characters are not equal:

        >>> msa2 = TabularMSA([DNA('ACG'), DNA('--G')])
        >>> msa1 != msa2
        True

        MSAs with different types of sequences (different ``dtype``) are not
        equal:

        >>> msa3 = TabularMSA([RNA('ACG'), RNA('AC-')])
        >>> msa1 != msa3
        True

        MSAs with different sequence metadata are not equal:

        >>> msa4 = TabularMSA([DNA('ACG', metadata={'id': 'a'}), DNA('AC-')])
        >>> msa1 != msa4
        True

        MSAs with different keys are not equal:

        >>> msa5 = TabularMSA([DNA('ACG'), DNA('AC-')], key=str)
        >>> msa1 != msa5
        True

        """
        return not (self == other)

    @experimental(as_of='0.4.0-dev')
    def has_keys(self):
        """Determine if keys exist on the MSA.

        Returns
        -------
        bool
            Indicates whether the MSA has keys.

        See Also
        --------
        keys
        reindex

        Examples
        --------
        >>> from skbio import DNA, TabularMSA
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')])
        >>> msa.has_keys()
        False
        >>> msa = TabularMSA([DNA('ACG'), DNA('AC-')], key=str)
        >>> msa.has_keys()
        True

        """
        return self._keys is not None

    @experimental(as_of='0.4.0-dev')
    def reindex(self, key=None, keys=None):
        """Reassign keys to sequences in the MSA.

        Parameters
        ----------
        key : callable or metadata key, optional
            If provided, defines a unique, hashable key for each sequence in
            the MSA. Can either be a callable accepting a single argument (each
            sequence) or a key into each sequence's ``metadata`` attribute.
        keys : iterable, optional
            An iterable of the same length as the number of sequences in the
            MSA. `keys` must contain unique, hashable elements. Each element
            will be used as the respective key for the sequences in the MSA.

        Raises
        ------
        ValueError
            If `key` and `keys` are both provided.
        ValueError
            If `keys` is not the same length as the number of sequences in the
            MSA.
        UniqueError
            If keys are not unique.

        See Also
        --------
        keys
        has_keys

        Notes
        -----
        If `key` or `keys` are not provided, keys will not be set and certain
        operations requiring keys will raise an ``OperationError``.

        Examples
        --------
        Create a ``TabularMSA`` object without keys:

        >>> from skbio import DNA, TabularMSA
        >>> seqs = [DNA('ACG', metadata={'id': 'a'}),
        ...         DNA('AC-', metadata={'id': 'b'})]
        >>> msa = TabularMSA(seqs)
        >>> msa.has_keys()
        False

        Set keys on the MSA, using each sequence's ID:

        >>> msa.reindex(key='id')
        >>> msa.has_keys()
        True
        >>> msa.keys # doctest: +NORMALIZE_WHITESPACE
        array(['a', 'b'],
              dtype='<U1')

        Remove keys from the MSA:

        >>> msa.reindex()
        >>> msa.has_keys()
        False

        Alternatively, an iterable of keys may be passed via `keys`:

        >>> msa.reindex(keys=['a', 'b'])
        >>> msa.keys # doctest: +NORMALIZE_WHITESPACE
        array(['a', 'b'],
              dtype='<U1')

        """
        if key is not None and keys is not None:
            raise ValueError(
                "Cannot use both `key` and `keys` at the same time.")

        keys_ = None
        if key is not None:
            keys_ = [resolve_key(seq, key) for seq in self._seqs]
        elif keys is not None:
            keys = list(keys)
            if len(keys) != len(self):
                raise ValueError(
                    "Number of elements in `keys` must match number of "
                    "sequences: %d != %d" % (len(keys), len(self)))
            keys_ = keys

        if keys_ is not None:
            # Hashability of keys is implicitly checked here.
            duplicates = find_duplicates(keys_)
            if duplicates:
                raise UniqueError(
                    "Keys must be unique. Duplicate keys: %r" % duplicates)
            # Create an immutable ndarray to ensure key invariants are
            # preserved.
            keys_ = np.asarray(keys_)
            keys_.flags.writeable = False

        self._keys = keys_

    @experimental(as_of='0.4.0-dev')
    def append(self, seq):
        if type(seq) != self._dtype:
            raise TypeError("Appended sequence type must be %s" %
                            self._dtype)
        elif len(seq) != self.shape.position:
            raise TypeError("Appended sequence length must be %d, but was %d" %
                            (self.shape.position, len(seq)))
        else:
            self._shape = _Shape(sequence=self._shape.sequence+1,
                                 position=self._shape.position)
            self._seqs.append(seq)
