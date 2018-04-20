# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


class IOSourceError(Exception):
    """Raised when a file source cannot be resolved."""
    pass


class FileFormatError(Exception):
    """Raised when a file cannot be parsed."""
    pass


class UnrecognizedFormatError(FileFormatError):
    """Raised when a file's format is unknown, ambiguous, or unidentifiable."""
    pass


class GenBankFormatError(FileFormatError):
    """Raised when a ``genbank`` formatted file cannot be parsed."""
    pass


class EMBLFormatError(FileFormatError):
    """Raised when a ``EMBL`` formatted file cannot be parsed."""
    pass


class GFF3FormatError(FileFormatError):
    """Raised when a ``GFF3`` formatted file cannot be parsed."""
    pass


class BLAST7FormatError(FileFormatError):
    """Raised when a ``blast7`` formatted file cannot be parsed."""
    pass


class ClustalFormatError(FileFormatError):
    """Raised when a ``clustal`` formatted file cannot be parsed."""
    pass


class FASTAFormatError(FileFormatError):
    """Raised when a ``fasta`` formatted file cannot be parsed."""
    pass


class QUALFormatError(FASTAFormatError):
    """Raised when a ``qual`` formatted file cannot be parsed."""
    pass


class LSMatFormatError(FileFormatError):
    """Raised when a ``lsmat`` formatted file cannot be parsed."""
    pass


class OrdinationFormatError(FileFormatError):
    """Raised when an ``ordination`` formatted file cannot be parsed."""
    pass


class NewickFormatError(FileFormatError):
    """Raised when a ``newick`` formatted file cannot be parsed."""
    pass


class FASTQFormatError(FileFormatError):
    """Raised when a ``fastq`` formatted file cannot be parsed."""
    pass


class PhylipFormatError(FileFormatError):
    """Raised when a ``phylip`` formatted file cannot be parsed.

    May also be raised when an object (e.g., ``TabularMSA``) cannot be written
    in ``phylip`` format.

    """
    pass


class QSeqFormatError(FileFormatError):
    """Raised when a ``qseq`` formatted file cannot be parsed."""
    pass


class StockholmFormatError(FileFormatError):
    """Raised when a ``stockholm`` formatted file cannot be parsed."""
    pass


class InvalidRegistrationError(Exception):
    """Raised if function doesn't meet the expected API of its registration."""
    pass


class DuplicateRegistrationError(Exception):
    """Raised when a function is already registered in skbio.io"""
    pass
