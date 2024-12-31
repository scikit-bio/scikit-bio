# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re


def _note_into_doc(note, doc):
    """Insert a note into a docstring under the head line.

    Parameters
    ----------
    note : str
        Note to insert.
    doc : str
        Docstring to modify.

    Returns
    -------
    str
        Modified docstring.

    Notes
    -----
    The note will be placed below the first paragraph (usually a single line) of the
    docstring, like a subtitle.

    """
    # Make sure note has one line end.
    note = note.rstrip() + "\n"

    # No docstring: note becomes the entire docstring.
    if not doc:
        return note

    # Determine the indentation of the second non-empty line.
    # Correct indentation is important for docstring rendering.
    if match := re.search(r"\n+(\s*)\S*", doc):
        indent = match.group(1)

    # If such a line doesn't exist, i.e., the docstring is a single line, then
    # indentation is empty.
    else:
        indent = ""

    # Find the first paragraph break (an empty line) and insert the note after the
    # first paragraph.
    if pos := doc.find("\n\n") + 1:
        return f"{doc[:pos].rstrip()}\n\n{indent}{note}{doc[pos:]}"

    # If not found, i.e., the docstring is a single paragraph, then append the note to
    # the end of docstring.
    else:
        return f"{doc.rstrip()}\n\n{indent}{note}{indent}"


def _note_into_doc_param(note, doc, param):
    """Insert a note into a docstring under the description of a parameter.

    Parameters
    ----------
    note : str
        Note to insert.
    doc : str
        Docstring to modify.
    param : str
        Target parameter.

    Returns
    -------
    str
        Modified docstring.

    Raises
    ------
    ValueError
        Parameter is missing or its format is invalid.

    Notes
    -----
    The note will be placed below the description of the parameter, and above the next
    parameter or the next section.

    """
    note = note.rstrip() + "\n"

    # Find the header line of the parameter.
    if not (match := re.search(rf"(^\s*){param}\s*:.*?\n", doc, re.MULTILINE)):
        raise ValueError(
            f'Parameter "{param}" is missing from the docstring or its format is '
            "invalid."
        )

    # Determine the indentation of the header line.
    indent = match.group(1)

    # Find the next line with the same or less indentation, which indicates the next
    # parameter or the end of the "Parameters" section.
    after = doc[(end := match.end()) :]
    match = re.search(rf"^{indent}(\S|$)", after, re.MULTILINE)

    # Determine the insertion point.
    pos = end + (match.start() if match else len(after))

    # Insert the message. It should have 1+ indentation level (i.e., 4 spaces) than
    # the header line.
    return f"{doc[:pos].rstrip()}\n\n{indent}    {note}\n{doc[pos:]}"


def _deprecation_note(ver=None, msg=None):
    """Create a note indicating deprecation."""
    note = f".. deprecated:: {ver}"
    if msg:
        note += " " + msg
    return note


def _renaming_note(name, ver=None, warn=False):
    """Create a message indicating renaming."""
    if not ver:
        return f"Alias: ``{name}``"
    note = (
        f".. versionchanged:: {ver} "
        f"Renamed from ``{name}``. The old name is kept as an alias"
    )
    if warn:
        note += " but is deprecated"
    note += "."
    return note
