# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# Script to manipulate sitemap.xml.

import os
import re


# -- Configuration -----------------------------------------------------------

rootdir = "build/html"
sitemap = 'sitemap.xml'


# -- Workflow ----------------------------------------------------------------

pattern = re.compile(r'<url>(.*?)</url>', flags=re.DOTALL)


def substitute(match):
    s = match.group()

    # omit pages like "__eq__.html"
    if '__.html' in s:
        return ''

    return s


cwd = os.getcwd()
os.chdir(os.path.join(os.path.dirname(__file__), rootdir))

with open(sitemap, 'r') as f:
    content = f.read()

content = pattern.sub(substitute, content)

with open(sitemap, 'w') as f:
    f.write(content)

os.chdir(cwd)
