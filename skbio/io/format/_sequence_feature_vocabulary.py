# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re

from skbio.metadata import IntervalMetadata
from skbio.io.format._base import _line_generator
from skbio.io import FileFormatError


def _vocabulary_change(format="insdc", read_in=True):
    """Return a dict that converts between memory and output vocabulary."""
    convert = {
        "phase": {"insdc": "codon_start"},
        "source": {"insdc": "inference"},
        "db_xref": {"gff3": "Dbxref"},
        "note": {"gff3": "Note"},
    }
    if read_in:
        return {v[format]: k for k, v in convert.items() if format in v}
    else:
        return {k: v[format] for k, v in convert.items() if format in v}


def _vocabulary_skip(format="insdc"):
    """Return vocabluaries skipped for auto disk output, given a format.

    Return a list of vocabularies that should be skipped when auto output to disk
    for the specified format.

    """
    skip = {
        "type": ("insdc", "gff3"),
        "ID": ("insdc"),
        "translation": ("gff3"),
        "strand": ("insdc"),
    }
    return [k for k, v in skip.items() if format in v]


def _yield_section(is_another_section, **kwargs):
    """Return function that returns successive sections from file.

    Parameters
    ----------
    is_another_section : callable
        It takes a string as input and return a boolean indicating
        a new section starts.
    kwargs : dict, optional
        Keyword arguments will be passed to `_line_generator`.

    Returns
    -------
    function
        A function accept a list of lines as input and return
        a generator to yield section one by one.

    """

    def parser(lines):
        curr = []
        for line in _line_generator(lines, **kwargs):
            # if we find another, return the previous section
            if is_another_section(line):
                if curr:
                    yield curr
                    curr = []
            curr.append(line)
        # don't forget to return the last section in the file
        if curr:
            yield curr

    return parser


def _parse_section_default(
    lines, label_delimiter=None, join_delimiter=" ", return_label=False
):
    """Parse sections in default way.

    Do 2 things:
        1. split first line with label_delimiter for label
        2. join all the lines into one str with join_delimiter.
    """
    data = []
    label = None
    line = lines[0]

    items = line.split(label_delimiter, 1)

    if len(items) == 2:
        label, section = items
    else:
        label = items[0]
        section = ""
    data.append(section)

    data.extend(lines[1:])
    data = join_delimiter.join(i.strip() for i in data)
    if return_label:
        return label, data
    else:
        return data


def _serialize_section_default(header, obj, indent=12):
    return "{header:<{indent}}{obj}\n".format(header=header, obj=obj, indent=indent)


def _parse_feature_table(lines, length):
    """Parse DDBJ/ENA/GenBank Feature Table."""
    imd = IntervalMetadata(length)
    # skip the 1st FEATURES line
    if lines[0].startswith("FEATURES"):
        lines = lines[1:]
    # magic number 21: the lines following header of each feature
    # are indented with 21 spaces.
    feature_indent = " " * 21
    section_splitter = _yield_section(
        lambda x: not x.startswith(feature_indent), skip_blanks=True, strip=False
    )

    for section in section_splitter(lines):
        _parse_single_feature(section, imd)
    return imd


def _parse_single_feature(lines, imd):
    """Parse a feature and add it to ``IntervalMetadata`` object.

    Parameters
    ----------
    imd : IntervalMetadata
        An IntervalMetadata object to which the parsed feature will be added.
    lines : list of strings
        A list of strings representing the lines of text to be parsed.

    """
    voca_change = _vocabulary_change("insdc")

    # each component of a feature starts with '/', except the 1st
    # component of location.
    section_splitter = _yield_section(lambda x: x.startswith("/"), strip=True)
    section_iter = section_splitter(lines)

    # 1st section is location
    section = next(section_iter)
    feature_type, feature_loc = _parse_section_default(
        section, join_delimiter="", return_label=True
    )

    metadata = {"type": feature_type, "__location": feature_loc}

    intvl = imd.add(*_parse_loc_str(feature_loc))

    for section in section_iter:
        # following sections are Qualifiers
        k, v = _parse_section_default(
            section, label_delimiter="=", join_delimiter=" ", return_label=True
        )
        # 1st char is '/'
        k = k[1:]
        if k in voca_change:
            k = voca_change[k]

        if k == "phase":
            v = int(v) - 1

        # some Qualifiers can appear multiple times
        if k in metadata:
            if not isinstance(metadata[k], list):
                metadata[k] = [metadata[k]]
            metadata[k].append(v)
        else:
            metadata[k] = v

    intvl.metadata.update(metadata)


def _parse_loc_str(loc_str):
    """Parse location string.

    .. warning: This converts coordinates to 0-based from 1-based
    GenBank coordinate system.

    The location descriptor can be one of the following [1]_:
    (a) a single base number. e.g. 467
    (b) a site between two indicated adjoining bases. e.g. 123^124
    (c) a single base chosen from within a specified range of bases (not
        allowed for new entries). e.g. 102.110
    (d) the base numbers delimiting a sequence span. e.g.340..565
    (e) a remote entry identifier followed by a local location
        descriptor (i.e., a-d). e.g. J00194.1:100..202

    Notes
    -----
    This does not fully handle (e) case. It will discard the remote
    entry part and only keep the local part. When it parses locations
    across strand (e.g. "complement(123..145),200..209"), it will
    record all the span parts but will record strand as negative.

    References
    ----------
    .. [1] http://www.insdc.org/files/feature_table.html#3.4

    """
    # define the tokens
    operators = ["join", "complement", "order"]
    LPAREN = r"(?P<LPAREN>\()"
    RPAREN = r"(?P<RPAREN>\))"
    COMMA = r"(?P<COMMA>,)"
    WS = r"(?P<WS>\s+)"
    a = r"(?P<A>\d+)"
    b = r"(?P<B>\d+\^\d+)"
    c = r"(?P<C>\d+\.\d+)"
    d = r"(?P<D><?\d+\.\.>?\d+)"
    e_left = r"(?P<EL><?[a-zA-Z_0-9\.]+:\d+\.\.>?\d+)"
    e_right = r"(?P<ER><?\d+\.\.>?[a-zA-Z_0-9\.]+:\d+)"
    illegal = r"(?P<ILLEGAL>.+)"
    # The order of tokens in the master regular expression also
    # matters. When matching, re tries to match pattens in the order
    # specified. Thus, if a pattern happens to be a substring of a
    # longer pattern, you need to make sure the longer pattern goes
    # first.
    master_pat = re.compile(
        "|".join(
            operators
            + [WS, LPAREN, RPAREN, COMMA, b, c, d, e_left, e_right, a, illegal]
        )
    )

    scanner = master_pat.scanner(loc_str)

    bounds = []
    fuzzy = []

    metadata = {"strand": "+"}

    for m in iter(scanner.match, None):
        p, v = m.lastgroup, m.group()
        if v == "complement":
            metadata["strand"] = "-"
        elif p == "A":
            start = int(v)
            bounds.append((start - 1, start))
            fuzzy.append((False, False))
        elif p == "B":
            start, end = v.split("^")
            start = int(start)
            bounds.append((start - 1, start))
            fuzzy.append((False, False))
        elif p == "C" or p == "D":
            if p == "C":
                start, end = v.split(".")
            else:
                start, end = v.split("..")
            fuzzy_s = fuzzy_e = False
            if start.startswith("<"):
                start = start[1:]
                fuzzy_s = True
            if end.startswith(">"):
                end = end[1:]
                fuzzy_e = True
            bounds.append((int(start) - 1, int(end)))
            fuzzy.append((fuzzy_s, fuzzy_e))
        elif p == "ILLEGAL":
            raise FileFormatError('Could not parse location string: "%s"' % loc_str)

    return bounds, fuzzy, metadata


def _serialize_feature_table(intervals, indent=21):
    """Serialize a list of intervals into a feature table format.

    Parameters
    ----------
    intervals : list of ``Interval``
        A list of Interval objects representing the intervals to be serialized.
    indent : int, optional
        The number of spaces to indent each serialized feature. Defaults to 21.

    """
    for intvl in intervals:
        yield _serialize_single_feature(intvl, indent)


def _serialize_single_feature(intvl, indent=21):
    """Serialize a single interval into feature format.

    Parameters
    ----------
    intvl : Interval
        The Interval object representing the interval to be serialized.
    indent : int, optional
        The number of spaces to indent each serialized feature. Defaults to 21.

    """
    # there are 5 spaces before Feature Key starts.
    padding = " " * 5
    qualifiers = []
    md = intvl.metadata
    voca_skip = _vocabulary_skip("insdc")
    voca_change = _vocabulary_change("insdc", read_in=False)
    # sort it so the output order is deterministic
    for k in sorted(md):
        if k.startswith("__") or k in voca_skip:
            continue
        v = md[k]
        if k == "phase":
            v = str(v + 1)
        if k in voca_change:
            k = voca_change[k]
        if isinstance(v, list):
            for vi in v:
                qualifiers.append(_serialize_qualifier(k, vi))
        else:
            qualifiers.append(_serialize_qualifier(k, v))

    if "__location" in md:
        loc = md["__location"]
    else:
        loc = _serialize_location(intvl)
    # the qualifiers start at column 22
    qualifiers = [" " * indent + i for i in qualifiers]
    return "{header:<{indent}}{loc}\n{qualifiers}\n".format(
        header=padding + md["type"],
        loc=loc,
        indent=indent,
        qualifiers="\n".join(qualifiers),
    )


def _serialize_location(intvl):
    loc = []
    for bound, fuzzy in zip(intvl.bounds, intvl.fuzzy):
        start, end = bound
        start += 1
        if start == end:
            s = str(start)
        elif fuzzy[0] and fuzzy[1]:
            s = "<%d..>%d" % (start, end)
        elif fuzzy[0] and not fuzzy[1]:
            s = "<%d..%d" % (start, end)
        elif not fuzzy[0] and fuzzy[1]:
            s = "%d..>%d" % (start, end)
        else:
            s = "%d..%d" % (start, end)
        loc.append(s)
    if len(loc) > 1:
        loc_str = "join({})".format(",".join(loc))
    else:
        loc_str = loc[0]
    if intvl.metadata.get("strand") == "-":
        loc_str = "complement({})".format(loc_str)
    return loc_str


def _serialize_qualifier(key, value):
    """Serialize a Qualifier in a feature.

    Parameters
    ----------
    key : str
        The key of the Qualifier, representing the type or name of the information.
    value : int, str
        The value associated with the Qualifier.

    """
    # if value is empty
    if not value:
        return "/%s" % key

    return "/{k}={v}".format(k=key, v=value)
