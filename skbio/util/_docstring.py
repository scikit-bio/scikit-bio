# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Docstring utilities."""

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


_BACKEND_DEVICES = {
    "numpy": ["cpu"],
    "cupy": ["gpu"],
    "torch": ["cpu", "gpu"],
    "jax": ["cpu", "gpu"],
    "dask": ["cpu"],
}

_BACKEND_DISPLAY = {
    "numpy": "NumPy",
    "cupy": "CuPy",
    "torch": "PyTorch",
    "jax": "JAX",
    "dask": "Dask",
}


def _array_api_compat_section(backends, devices=None):
    """Generate an RST compatibility table for the Array API notes section.

    Parameters
    ----------
    backends : list of str
        Supported backends.
    devices : list of str, optional
        Supported devices. Defaults to per-backend defaults.

    Returns
    -------
    str
        RST-formatted compatibility table paragraph.

    """
    # Collect the set of all devices to display as columns
    if devices is not None:
        all_devices = list(dict.fromkeys(devices))
    else:
        all_devices = []
        for b in backends:
            for d in _BACKEND_DEVICES.get(b, ["cpu"]):
                if d not in all_devices:
                    all_devices.append(d)

    col_headers = [d.upper() for d in all_devices]

    # Build rows
    rows = []
    for b in backends:
        display = _BACKEND_DISPLAY.get(b, b.capitalize())
        supported = _BACKEND_DEVICES.get(b, ["cpu"])
        if devices is not None:
            supported = [d for d in devices if d in supported]
        cells = ["\u2713" if d in supported else "n/a" for d in all_devices]
        rows.append((display, cells))

    # Compute column widths
    backend_col_w = max(len("Backend"), max(len(r[0]) for r in rows))
    device_col_ws = [
        max(len(h), max(len(r[1][i]) for r in rows)) for i, h in enumerate(col_headers)
    ]

    def hline(left, mid, right, fill="-"):
        parts = [fill * (backend_col_w + 2)]
        for w in device_col_ws:
            parts.append(fill * (w + 2))
        return "+" + "+".join(parts) + "+"

    def row_line(backend, cells, fill=" "):
        parts = [f" {backend:<{backend_col_w}} "]
        for val, w in zip(cells, device_col_ws):
            parts.append(f" {val:<{w}} ")
        return "|" + "|".join(parts) + "|"

    header_sep = hline("+", "+", "+", "=")
    normal_sep = hline("+", "+", "+", "-")
    header_row = row_line("Backend", col_headers)

    lines = [normal_sep, header_row, header_sep]
    for display, cells in rows:
        lines.append(row_line(display, cells))
        lines.append(normal_sep)
    table = "\n".join(lines)

    intro = (
        "This function supports the `Python array API standard"
        " <https://data-apis.org/array-api/latest/>`_."
        " Compatible array backends:\n\n"
    )
    return intro + table + "\n"


def _insert_into_notes_section(note, doc):
    """Insert a note into the Notes section of a docstring.

    Parameters
    ----------
    note : str
        Text to insert.
    doc : str
        Docstring to modify.

    Returns
    -------
    str
        Modified docstring with note prepended inside Notes section.

    Notes
    -----
    If a Notes section is found, the note is prepended inside it. Otherwise,
    a new Notes section is appended at the end of the docstring.

    """
    note = note.rstrip() + "\n"

    # Determine indentation from the second non-empty line (same as _note_into_doc)
    if match := re.search(r"\n+(\s*)\S*", doc):
        indent = match.group(1)
    else:
        indent = ""

    # Look for an existing Notes section "Notes\n-----"
    notes_match = re.search(r"^(\s*)Notes\s*\n\s*-{3,}", doc, re.MULTILINE)
    if notes_match:
        # Insert after the section header (and any blank lines following it)
        header_end = notes_match.end()
        # Skip any leading blank lines inside the section
        rest = doc[header_end:]
        skip = re.match(r"\n*", rest)
        insert_pos = header_end + (skip.end() if skip else 0)
        # Re-indent the note to match section body indentation
        indented_note = "".join(
            indent + line if line.strip() else line
            for line in note.splitlines(keepends=True)
        )
        return doc[:insert_pos] + indented_note + "\n" + doc[insert_pos:]

    # No Notes section: append one at the end
    doc_stripped = doc.rstrip()
    notes_section = f"\n\n{indent}Notes\n{indent}-----\n{indent}" + note.rstrip() + "\n"
    return doc_stripped + notes_section


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
