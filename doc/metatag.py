# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Script to insert a meta description tag into each documentation page. The
# content is the first line of the docstring of the corresponding object.

import os
import re
import glob
import html


# -- Configuration -----------------------------------------------------------

rootdir = "build/html"
max_length = 160


# -- Workflow ----------------------------------------------------------------

dd_pattern = re.compile(r'<dd>(.*?)</dd>', flags=re.DOTALL)
p_pattern = re.compile(r'<p>(.*?)</p>', flags=re.DOTALL)
meta_pattern = re.compile(r'(\s*)(<meta\b[^>]*\/?>)')

cwd = os.getcwd()
os.chdir(os.path.join(os.path.dirname(__file__), rootdir))

for file in glob.glob("**/*.html", recursive=True):
    with open(file, "r") as fh:
        content = fh.read()

    # find first line of docstring (summary)
    dd_match = dd_pattern.search(content)
    if not dd_match:
        continue
    p_match = p_pattern.search(dd_match.group(1))
    if not p_match:
        continue
    summary = p_match.group(1)

    # truncate summary to given length
    if len(summary) > max_length:
        summary = summary[:max_length - 3]
        if ' ' in summary:
            summary = summary.rsplit(' ', 1)[0]
        summary += "..."

    # make summary safe for HTML
    summary = html.escape(summary)

    # insert summary into metadata
    line = f'<meta name="description" content="{summary}" />'
    meta_match = meta_pattern.search(content)
    if not meta_match:
        continue
    indent = meta_match.group(1)
    pos = meta_match.end()
    content = content[:pos] + indent + line + content[pos:]

    with open(file, "w") as fh:
        fh.write(content)

os.chdir(cwd)
