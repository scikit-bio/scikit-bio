This README file describes the FASTQ example files provided as supplementary
information to the open-access publication:

P.J.A. Cock, C.J. Fields, N. Goto, M.L. Heuer and P.M. Rice (2009). The Sanger
FASTQ file format for sequences with quality scores, and the Solexa/Illumina
FASTQ variants.

These files are provided freely and we encourage anyone writing a FASTQ parser
to use them as part of your test suite. Permission is granted to freely
distribute and modify the files. We request (but do not insist) that this
README file is included, or at least a reference to the above paper. Please
cite the above paper if appropriate. We also request (but do not insist) that
the example files are not modified, in order that they may serve as a common
reference.

Invalid FASTQ files
===================

The archive contains the following sample FASTQ files with names of the form
error_NAME.fastq, which all contain errors and should be rejected (if parsed
as any of the three FASTQ variants):

error_diff_ids.fastq
error_double_qual.fastq
error_double_seq.fastq
error_long_qual.fastq
error_no_qual.fastq
error_qual_del.fastq
error_qual_escape.fastq
error_qual_null.fastq
error_qual_space.fastq
error_qual_tab.fastq
error_qual_unit_sep.fastq
error_qual_vtab.fastq
error_short_qual.fastq
error_spaces.fastq
error_tabs.fastq
error_trunc_at_seq.fastq
error_trunc_at_plus.fastq
error_trunc_at_qual.fastq
error_trunc_in_title.fastq
error_trunc_in_seq.fastq
error_trunc_in_plus.fastq
error_trunc_in_qual.fastq

Of these, those with names error_qual_XXX.fastq would be valid except for the
inclusion of spaces or non-printing ASCII characters outside the range allowed
in the quality string. The files named error_trunc_XXX.fastq would be valid
but for being truncated (e.g. simulating a partial copy over the network).

The special cases of FASTQ files which would be valid as one variant, but not
another, are covered below.

Valid FASTQ
===========

The archive contains the following valid sample FASTQ input files for testing:

longreads_original_sanger.fastq
wrapping_original_sanger.fastq
illumina_full_range_original_illumina.fastq
sanger_full_range_original_sanger.fastq
solexa_full_range_original_solexa.fastq
misc_dna_original_sanger.fastq
misc_rna_original_sanger.fastq

These all have the form NAME_original_FORMAT.fastq, where NAME is a prefix for
that example, and FORMAT is one of sanger, solexa or illumina indicating which
FASTQ variant that example is using. There are three matching files called
NAME_as_FORMAT.fastq showing how the original file should be converted into
each of the three FASTQ variants. These converted files are standardised not
to use line wrapping (so each record has exactly four lines), and omit the
optional repetition of the read titles on the plus line.

The file longreads_original_sanger.fastq is based on real Roche 454 reads from
the Sanger Institute for the the potato cyst nematodes Globodera pallida. Ten
of the reads have been presented as FASTQ records, wrapping the sequence and
the quality lines at 80 characters. This means some of the quality lines start
with "@" or "+" characters, which may cause problems with naive parsers. Also
note that the sequence is mixed case (with upper case denoting the trimmed
region), and furthermore the free format title lines are over 100 characters
and encode assorted read information (and are repeated on the "+" lines).

The wrapping_original_sanger.fastq is based on three real reads from the NCBI
Short Read Archive, but has been carefully edited to use line wrapping for the
quality lines (but not the sequence lines) such that the due to the occurrence
of "@" and "+" on alternating lines, the file may be misinterpreted by a
simplistic parser. While this is therefore a very artificial example, it
remains a valid FASTQ file, and is useful for testing purposes.

The sanger_full_range_original_sanger.fastq file uses PHRED scores from 0 to
93 inclusive, covering ASCII characters from 33 (!) to 126 (~). This means it
cannot be treated as a Solexa or Illumina 1.3+ FASTQ file, and attempting to
parse it as such should raise an error.

The solexa_full_range_original_solexa.fastq file uses Solexa scores from -5 to
62 inclusive, covering ASCII characters from 59 (;) to 126 (~). This means it
cannot be treated as a Illumina 1.3+ FASTQ file, and attempting to parse it as
such should raise an error. On the basis of the quality characters, the file
would also qualify as a valid Sanger FASTQ file.

The illumina_full_range_original_illumina.fastq file uses PHRED scores from 0
to 62 inclusive, covering ASCII characters from 64 (@) to 126 (~). On the
basis of the quality characters, the file would also qualify as a valid Sanger
or Solexa FASTQ file.

The misc_dna_original_sanger.fastq and misc_rna_original_sanger.fastq files
are artificial reads using the full range of IUPAC DNA or RNA letters,
including ambiguous character codes, and both cases.
