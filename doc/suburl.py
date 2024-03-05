# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Script to substitute domain names with relative paths in hyperlinks in the
# built html files, such as to avoid unnecessary domain resolutions when the
# visitor navigates the website. It may be executed only when "doc" and "web"
# are hosted at the same web server.

import os
import re
import glob
from functools import partial


# -- Configuration -----------------------------------------------------------

rootdir = "build/html"
source = "https://scikit.bio"
target = "../.."


# -- Workflow ----------------------------------------------------------------

pattern = re.compile(f'href="{re.escape(source)}/([^"]+)"')


def substitute(match, prefix):
    return f'href="{prefix}{target}/{match.group(1)}"'


cwd = os.getcwd()
os.chdir(os.path.join(os.path.dirname(__file__), rootdir))

for file in glob.glob("**/*.html", recursive=True):
    depth = len(os.path.normpath(file).split(os.sep))
    prefix = "../" * (depth - 1)
    with open(file, "r") as fh:
        content = fh.read()
    content = content.replace(
        f'href="{source}"', f'href="{prefix}{target}/index.html"'
    )
    repl = partial(substitute, prefix=prefix)
    content = pattern.sub(repl, content)
    with open(file, "w") as fh:
        fh.write(content)

os.chdir(cwd)
