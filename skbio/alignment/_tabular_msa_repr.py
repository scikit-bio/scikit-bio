# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import six

import itertools
import numbers
import textwrap

from skbio.sequence._base import ElasticLines


class _TabularMSAReprBuilder(object):
    """Build a ``TabularMSA`` repr.

    Parameters
    ----------
    seq : Sequence
        Sequence to repr.
    width : int
        Maximum width of the repr.
    indent : int
        Number of spaces to use for indented lines.
    chunk_size: int
        Number of characters in each chunk of a sequence.

    """
    def __init__(self, seq, width, indent, chunk_size):
        self._msa = seq
        self._width = width
        self._indent = ' ' * indent
        self._chunk_size = chunk_size

    def build(self):
        lines = ElasticLines()

        cls_name = self._msa.__class__.__name__
        dtype_class_name = self._msa._seqs[0].__class__.__name__
        lines.add_line(cls_name + ' <' + dtype_class_name + '>')
        lines.add_separator()

        if self._msa.has_metadata():
            lines.add_line('Metadata:')
            # Python 3 doesn't allow sorting of mixed types so we can't just
            # use sorted() on the metadata keys. Sort first by type then sort
            # by value within each type.
            for key in self._sorted_keys_grouped_by_type(self._msa.metadata):
                value = self._msa.metadata[key]
                lines.add_lines(self._format_metadata_key_value(key, value))

        if self._msa.has_positional_metadata():
            lines.add_line('Positional metadata:')
            for key in self._msa.positional_metadata.columns.values.tolist():
                dtype = self._msa.positional_metadata[key].dtype
                lines.add_lines(
                    self._format_positional_metadata_column(key, dtype))

        lines.add_line('Stats:')
        for label, value in self._msa._repr_stats():
            lines.add_line('%s%s: %s' % (self._indent, label, value))
        lines.add_separator()

        num_lines = self._msa.shape.sequence

        # display entire sequence if we can, else display the first two and
        # last two lines separated by ellipsis
        if num_lines <= 5:
            lines.add_lines(self._format_lines(range(num_lines)))
        else:
            lines.add_lines(self._format_lines(range(2)))
            lines.add_line('...')
            lines.add_lines(self._format_lines(range(num_lines - 2,
                                                     num_lines)))

        return lines.to_str()

    def _format_lines(self, line_indices):
        num_chars_at_end = 4
        lines = []
        for line_index in line_indices:
            seq_str = str(self._msa._seqs[line_index])
            formatted_seq = (seq_str[0:-num_chars_at_end] + ' ... ' +
                             seq_str[-num_chars_at_end:])
            lines.append(formatted_seq)
        return lines

    def _sorted_keys_grouped_by_type(self, dict_):
        """Group keys within a dict by their type and sort within type."""
        type_sorted = sorted(dict_, key=self._type_sort_key)
        type_and_value_sorted = []
        for _, group in itertools.groupby(type_sorted, self._type_sort_key):
            type_and_value_sorted.extend(sorted(group))
        return type_and_value_sorted

    def _type_sort_key(self, key):
        return repr(type(key))

    def _format_metadata_key_value(self, key, value):
        """Format metadata key:value, wrapping across lines if necessary."""
        key_fmt = self._format_key(key)

        supported_type = True
        if isinstance(value, (six.text_type, six.binary_type)):
            # for stringy values, there may be u'' or b'' depending on the type
            # of `value` and version of Python. find the starting quote
            # character so that wrapped text will line up with that instead of
            # the string literal prefix character. for example:
            #
            #     'foo': u'abc def ghi
            #              jkl mno'
            value_repr = repr(value)
            extra_indent = 1
            if not (value_repr.startswith("'") or value_repr.startswith('"')):
                extra_indent = 2
        # handles any number, this includes bool
        elif value is None or isinstance(value, numbers.Number):
            value_repr = repr(value)
            extra_indent = 0
        else:
            supported_type = False

        if not supported_type or len(value_repr) > 140:
            value_repr = str(type(value))
            # extra indent of 1 so that wrapped text lines up past the bracket:
            #
            #     'foo': <type
            #             'dict'>
            extra_indent = 1

        return self._wrap_text_with_indent(value_repr, key_fmt, extra_indent)

    def _format_key(self, key):
        """Format metadata key.

        Includes initial indent and trailing colon and space:

            <indent>'foo':<space>

        """
        key_fmt = self._indent + repr(key)
        supported_types = (six.text_type, six.binary_type, numbers.Number,
                           type(None))
        if len(key_fmt) > (self._width / 2) or not isinstance(key,
                                                              supported_types):
            key_fmt = self._indent + str(type(key))
        return '%s: ' % key_fmt

    def _wrap_text_with_indent(self, text, initial_text, extra_indent):
        """Wrap text across lines with an initial indentation.

        For example:

            'foo': 'abc def
                    ghi jkl
                    mno pqr'

        <indent>'foo':<space> is `initial_text`. `extra_indent` is 1. Wrapped
        lines are indented such that they line up with the start of the
        previous line of wrapped text.

        """
        return textwrap.wrap(
            text, width=self._width, expand_tabs=False,
            initial_indent=initial_text,
            subsequent_indent=' ' * (len(initial_text) + extra_indent))

    def _format_positional_metadata_column(self, key, dtype):
        key_fmt = self._format_key(key)
        dtype_fmt = '<dtype: %s>' % str(dtype)
        return self._wrap_text_with_indent(dtype_fmt, key_fmt, 1)
