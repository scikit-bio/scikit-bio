#! /usr/bin/env python
from setuptools import setup

from ordination import __version__

setup(
    name='Ordination',
    version=__version__,
    packages=['ordination', 'ordination.test'],
    package_data={
        'ordination.test': ['data/*']
        },
    license='BSD',
    install_requires=['numpy', 'matplotlib'],
)
