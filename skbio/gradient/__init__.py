r"""Gradient analyses (:mod:`skbio.gradient`)
===============================================

.. currentmodule:: skbio.gradient

This module provides functionality for performing gradient analyses.
The algorithms included in this module mainly allows performing analysis of
volatility on time series data, but they can be applied to any data that
contains a gradient.

Classes
-------

.. autosummary::
   :toctree: generated/

   GradientANOVA
   AverageGradientANOVA
   TrajectoryGradientANOVA
   FirstDifferenceGradientANOVA
   WindowDifferenceGradientANOVA
   GroupResults
   CategoryResults
   GradientANOVAResults

Examples
--------
Assume we have the following coordinates:

>>> import numpy as np
>>> import pandas as pd
>>> from skbio.gradient import AverageGradientANOVA
>>> coord_data = {'PC.354': np.array([0.2761, -0.0341, 0.0633, 0.1004]),
...               'PC.355': np.array([0.2364, 0.2186, -0.0301, -0.0225]),
...               'PC.356': np.array([0.2208, 0.0874, -0.3519, -0.0031]),
...               'PC.607': np.array([-0.1055, -0.4140, -0.15, -0.116]),
...               'PC.634': np.array([-0.3716, 0.1154, 0.0721, 0.0898])}
>>> coords = pd.DataFrame.from_dict(coord_data, orient='index')

the following metadata map:

>>> metadata_map = {'PC.354': {'Treatment': 'Control', 'Weight': '60'},
...            'PC.355': {'Treatment': 'Control', 'Weight': '55'},
...            'PC.356': {'Treatment': 'Control', 'Weight': '50'},
...            'PC.607': {'Treatment': 'Fast', 'Weight': '65'},
...            'PC.634': {'Treatment': 'Fast', 'Weight': '68'}}
>>> metadata_map = pd.DataFrame.from_dict(metadata_map, orient='index')

and the following array with the proportion explained of each coord:

>>> prop_expl = np.array([25.6216, 15.7715, 14.1215, 11.6913, 9.8304])

Then to compute the average trajectory of this data:

>>> av = AverageGradientANOVA(coords, prop_expl, metadata_map,
...                     trajectory_categories=['Treatment'],
...                     sort_category='Weight')
>>> trajectory_results = av.get_trajectories()

Check the algorithm used to compute the trajectory_results:

>>> print(trajectory_results.algorithm)
avg

Check if we weighted the data or not:

>>> print(trajectory_results.weighted)
False

Check the results of one of the categories:

>>> print(trajectory_results.categories[0].category)
Treatment
>>> print(trajectory_results.categories[0].probability)
0.0118478282382

Check the results of one group of one of the categories:

>>> print(trajectory_results.categories[0].groups[0].name)
Control
>>> print(trajectory_results.categories[0].groups[0].trajectory)
[ 3.52199973  2.29597001  3.20309816]
>>> print(trajectory_results.categories[0].groups[0].info)
{'avg': 3.007022633956606}


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._gradient import (
    GradientANOVA,
    AverageGradientANOVA,
    TrajectoryGradientANOVA,
    FirstDifferenceGradientANOVA,
    WindowDifferenceGradientANOVA,
    GroupResults,
    CategoryResults,
    GradientANOVAResults,
)


__all__ = [
    "GradientANOVA",
    "AverageGradientANOVA",
    "TrajectoryGradientANOVA",
    "FirstDifferenceGradientANOVA",
    "WindowDifferenceGradientANOVA",
    "GroupResults",
    "CategoryResults",
    "GradientANOVAResults",
]
