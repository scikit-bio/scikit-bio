r"""
Gradient analyses (:mod:`skbio.stats.gradient`)
===============================================

.. currentmodule:: skbio.stats.gradient

This module provides functionality for performing gradient analyses.
The algorithms included in this module mainly allows performing analysis of
volatility on time series data, but they can be applied to any data that
contains a gradient.

Classes
-------

.. autosummary::
   :toctree:

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
>>> from skbio.stats.gradient import AverageGradientANOVA
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
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from copy import deepcopy
from collections import defaultdict
from numbers import Integral

import numpy as np
from natsort import realsorted
from scipy.stats import f_oneway

from skbio.util._decorator import experimental


def _weight_by_vector(trajectories, w_vector):
    r"""weights the values of `trajectories` given a weighting vector
    `w_vector`.

    Each value in `trajectories` will be weighted by the 'rate of change'
    to 'optimal rate of change' ratio. The 'rate of change' of a vector
    measures how each point in the vector changes with respect to its
    predecessor point. The 'optimal rate of change' is the rate of change
    in which each point in the vector performs the same change than its
    predecessor, meaning that when calling this function over evenly spaced
    `w_vector` values, no change will be reflected on the output.

    Parameters
    ----------
    trajectories: pandas.DataFrame
        Values to weight
    w_vector: pandas.Series
        Values used to weight `trajectories`

    Returns
    -------
    pandas.DataFrame
        A weighted version of `trajectories`.

    Raises
    ------
    ValueError
        If `trajectories` and `w_vector` don't have equal lengths
        If `w_vector` is not a gradient
    TypeError
        If `trajectories` and `w_vector` are not iterables
    """
    try:
        if len(trajectories) != len(w_vector):
            raise ValueError("trajectories (%d) & w_vector (%d) must be equal "
                             "lengths" % (len(trajectories), len(w_vector)))
    except TypeError:
        raise TypeError("trajectories and w_vector must be iterables")

    # check no repeated values are passed in the weighting vector
    if len(set(w_vector)) != len(w_vector):
        raise ValueError("The weighting vector must be a gradient")

    # no need to weight in case of a one element vector
    if len(w_vector) == 1:
        return trajectories

    # Cast to float so divisions have a floating point resolution
    total_length = float(max(w_vector) - min(w_vector))

    # Reflects the expected gradient between subsequent values in w_vector
    # the first value isn't weighted so subtract one from the number of
    # elements
    optimal_gradient = total_length/(len(w_vector)-1)

    # for all elements apply the weighting function
    for i, idx in enumerate(trajectories.index):
        # Skipping the first element is it doesn't need to be weighted
        if i != 0:
            trajectories.loc[idx] = (
                trajectories.loc[idx] * optimal_gradient /
                np.abs((w_vector[i] - w_vector[i-1]))
            )

    return trajectories


def _ANOVA_trajectories(category, res_by_group):
    r"""Run ANOVA over `res_by_group`

    If ANOVA cannot be run in the current category (because either there is
    only one group in category or there is a group with only one member)
    the result CategoryResults instance has `probability` and `groups` set
    to None and message is set to a string explaining why ANOVA was not run

    Returns
    -------
    CategoryResults
        An instance of CategoryResults holding the results of the trajectory
        analysis applied on `category`
    """
    # If there is only one group under category we cannot run ANOVA
    if len(res_by_group) == 1:
        return CategoryResults(category, None, None,
                               'Only one value in the group.')
    # Check if groups can be tested using ANOVA. ANOVA testing requires
    # all elements to have at least size greater to one.
    values = [res.trajectory.astype(float) for res in res_by_group]
    if any([len(value) == 1 for value in values]):
        return CategoryResults(category, None, None,
                               'This group can not be used. All groups '
                               'should have more than 1 element.')
    # We are ok to run ANOVA
    _, p_val = f_oneway(*values)
    return CategoryResults(category, p_val, res_by_group, None)


class GroupResults:
    """Store the trajectory results of a group of a metadata category

    Attributes
    ----------
    name : str
        The name of the group within the metadata category
    trajectory : array like
        The result trajectory in an 1-D numpy array
    mean : float
        The mean of the trajectory
    info : dict
        Any extra information computed by the trajectory algorithm. Depends on
        the algorithm
    message : str
        A message with information of the execution of the algorithm

    """

    @experimental(as_of="0.4.0")
    def __init__(self, name, trajectory, mean, info, message):
        self.name = name
        self.trajectory = trajectory
        self.mean = mean
        self.info = info
        self.message = message

    @experimental(as_of="0.4.0")
    def to_files(self, out_f, raw_f):
        r"""Save the trajectory analysis results for a category group to files
        in text format.

        Parameters
        ----------
        out_f : file-like object
            File-like object to write trajectory analysis data to. Must have a
            `write` method. It is the caller's responsibility to close
            `out_f` when done (if necessary)
        raw_f : file-like object
            File-like object to write trajectories trajectory values. Must have
            a `write` method. It is the caller's responsibility to close
            `out_f` when done (if necessary)
        """
        out_f.write('For group "%s", the group means is: %f\n'
                    % (self.name, self.mean))
        raw_f.write('For group "%s":\n' % self.name)

        if self.message:
            out_f.write('%s\n' % self.message)
            raw_f.write('%s\n' % self.message)

        out_f.write('The info is: %s\n'
                    % sorted(((k, v) for k, v in self.info.items())))
        raw_f.write('The trajectory is:\n[%s]\n'
                    % ", ".join(map(str, self.trajectory)))


class CategoryResults:
    """Store the trajectory results of a metadata category

    Attributes
    ----------
    category : str
        The name of the category
    probability : float
        The ANOVA probability that the category groups are independent
    groups : list of GroupResults
        The trajectory results for each group in the category
    message : str
        A message with information of the execution of the algorithm

    """

    @experimental(as_of="0.4.0")
    def __init__(self, category, probability, groups, message):
        self.category = category
        self.probability = probability
        self.groups = groups
        self.message = message

    @experimental(as_of="0.4.0")
    def to_files(self, out_f, raw_f):
        r"""Save the trajectory analysis results for a category to files in
        text format.

        Parameters
        ----------
        out_f : file-like object
            File-like object to write trajectory analysis data to. Must have a
            `write` method. It is the caller's responsibility to close `out_f`
            when done (if necessary)
        raw_f : file-like object
            File-like object to write trajectory raw values. Must have a
            `write` method. It is the caller's responsibility to close `out_f`
            when done (if necessary)
        """
        if self.probability is None:
            out_f.write('Grouped by "%s": %s\n'
                        % (self.category, self.message))
        else:
            out_f.write('Grouped by "%s", probability: %f\n'
                        % (self.category, self.probability))
            raw_f.write('Grouped by "%s"\n' % self.category)
            for group in self.groups:
                group.to_files(out_f, raw_f)


class GradientANOVAResults:
    """Store the trajectory results

    Attributes
    ----------
    algorithm : str
        The algorithm used to compute trajectories
    weighted : bool
        If true, a weighting vector was used
    categories : list of CategoryResults
        The trajectory results for each metadata category

    """

    @experimental(as_of="0.4.0")
    def __init__(self, algorithm, weighted, categories):
        self.algorithm = algorithm
        self.weighted = weighted
        self.categories = categories

    @experimental(as_of="0.4.0")
    def to_files(self, out_f, raw_f):
        r"""Save the trajectory analysis results to files in text format.

        Parameters
        ----------
        out_f : file-like object
            File-like object to write trajectories analysis data to. Must have
            a `write` method. It is the caller's responsibility to close
            `out_f` when done (if necessary)
        raw_f : file-like object
            File-like object to write trajectories raw values. Must have a
            `write` method. It is the caller's responsibility to close `out_f`
            when done (if necessary)
        """
        out_f.write('Trajectory algorithm: %s\n' % self.algorithm)
        raw_f.write('Trajectory algorithm: %s\n' % self.algorithm)

        if self.weighted:
            out_f.write('** This output is weighted **\n')
            raw_f.write('** This output is weighted **\n')

        out_f.write('\n')
        raw_f.write('\n')

        for cat_results in self.categories:
            cat_results.to_files(out_f, raw_f)
            out_f.write('\n')
            raw_f.write('\n')


class GradientANOVA:
    r"""Base class for the Trajectory algorithms

    Parameters
    ----------
    coords : pandas.DataFrame
        The coordinates for each sample id
    prop_expl : array like
        The numpy 1-D array with the proportion explained by each axis in
        coords
    metadata_map : pandas.DataFrame
        The metadata map, indexed by sample ids and columns are metadata
        categories
    trajectory_categories : list of str, optional
        A list of metadata categories to use to create the trajectories. If
        None is passed, the trajectories for all metadata categories are
        computed. Default: None, compute all of them
    sort_category : str, optional
        The metadata category to use to sort the trajectories. Default: None
    axes : int, optional
        The number of axes to account while doing the trajectory specific
        calculations. Pass 0 to compute all of them. Default: 3
    weighted : bool, optional
        If true, the output is weighted by the space between samples in the
        `sort_category` column

    Raises
    ------
    ValueError
        If any category of `trajectory_categories` is not present in
        `metadata_map`
        If `sort_category` is not present in `metadata_map`
        If `axes` is not between 0 and the maximum number of axes available
        If `weighted` is True and no `sort_category` is provided
        If `weighted` is True and the values under `sort_category` are not
        numerical
        If `coords` and `metadata_map` does not have samples in common
    """
    # Should be defined by the derived classes
    _alg_name = None

    @experimental(as_of="0.4.0")
    def __init__(self, coords, prop_expl, metadata_map,
                 trajectory_categories=None, sort_category=None, axes=3,
                 weighted=False):
        if not trajectory_categories:
            # If trajectory_categories is not provided, use all the categories
            # present in the metadata map
            trajectory_categories = metadata_map.keys()
        else:
            # Check that trajectory_categories are in metadata_map
            for category in trajectory_categories:
                if category not in metadata_map:
                    raise ValueError("Category %s not present in metadata."
                                     % category)

        # Check that sort_categories is in metadata_map
        if sort_category and sort_category not in metadata_map:
            raise ValueError("Sort category %s not present in metadata."
                             % sort_category)

        if axes == 0:
            # If axes == 0, we should compute the trajectories for all axes
            axes = len(prop_expl)
        elif axes > len(prop_expl) or axes < 0:
            # Axes should be 0 <= axes <= len(prop_expl)
            raise ValueError("axes should be between 0 and the max number of "
                             "axes available (%d), found: %d "
                             % (len(prop_expl), axes))

        # Restrict coordinates to those axes that we actually need to compute
        self._coords = coords.loc[:, :axes-1]
        self._prop_expl = prop_expl[:axes]
        self._metadata_map = metadata_map
        self._weighted = weighted

        # Remove any samples from coords not present in mapping file
        # and remove any samples from metadata_map not present in coords
        self._normalize_samples()

        # Create groups
        self._make_groups(trajectory_categories, sort_category)

        # Compute the weighting_vector
        self._weighting_vector = None
        if weighted:
            if not sort_category:
                raise ValueError("You should provide a sort category if you "
                                 "want to weight the trajectories")
            try:
                self._weighting_vector = \
                    self._metadata_map[sort_category].astype(np.float64)
            except ValueError:
                raise ValueError("The sorting category must be numeric")

        # Initialize the message buffer
        self._message_buffer = []

    @experimental(as_of="0.4.0")
    def get_trajectories(self):
        r"""Compute the trajectories for each group in each category and run
        ANOVA over the results to test group independence.

        Returns
        -------
        GradientANOVAResults
            An instance of GradientANOVAResults holding the results.
        """
        result = GradientANOVAResults(self._alg_name, self._weighted, [])
        # Loop through all the categories that we should compute
        # the trajectories
        for cat, cat_groups in self._groups.items():
            # Loop through all the category values present in the current
            # category and compute the trajectory for each of them
            res_by_group = []
            for group in sorted(cat_groups, key=lambda k: str(k)):
                res_by_group.append(
                    self._get_group_trajectories(group, cat_groups[group]))

            result.categories.append(_ANOVA_trajectories(cat, res_by_group))

        return result

    def _normalize_samples(self):
        r"""Ensures that `self._coords` and `self._metadata_map` have the same
        sample ids

        Raises
        ------
        ValueError
            If `coords` and `metadata_map` does not have samples in common
        """
        # Figure out the sample ids in common
        coords_sample_ids = set(self._coords.index)
        mm_sample_ids = set(self._metadata_map.index)
        sample_ids = coords_sample_ids.intersection(mm_sample_ids)

        # Check if they actually have sample ids in common
        if not sample_ids:
            raise ValueError("Coordinates and metadata map had no samples "
                             "in common")

        # Need to take a subset of coords
        if coords_sample_ids != sample_ids:
            self._coords = self._coords.loc[sample_ids]
        # Need to take a subset of metadata_map
        if mm_sample_ids != sample_ids:
            self._metadata_map = self._metadata_map.loc[sample_ids]

    def _make_groups(self, trajectory_categories, sort_category):
        r"""Groups the sample ids in `self._metadata_map` by the values in
        `trajectory_categories`

        Creates `self._groups`, a dictionary keyed by category and values are
        dictionaries in which the keys represent the group name within the
        category and values are ordered lists of sample ids

        If `sort_category` is not None, the sample ids are sorted based on the
        values under this category in the metadata map. Otherwise, they are
        sorted using the sample id.

        Parameters
        ----------
        trajectory_categories : list of str
            A list of metadata categories to use to create the groups.
            Default: None, compute all of them
        sort_category : str or None
            The category from self._metadata_map to use to sort groups
        """
        # If sort_category is provided, we used the value of such category to
        # sort. Otherwise, we use the sample id.
        if sort_category:
            def sort_val(sid):
                return self._metadata_map[sort_category][sid]
        else:
            def sort_val(sid):
                return sid

        self._groups = defaultdict(dict)
        for cat in trajectory_categories:
            # Group samples by category
            gb = self._metadata_map.groupby(cat)
            for g, df in gb:
                self._groups[cat][g] = realsorted(df.index, key=sort_val)

    def _get_group_trajectories(self, group_name, sids):
        r"""Compute the trajectory results for `group_name` containing the
        samples `sids`.

        Weights the data if `self._weighted` is True and ``len(sids) > 1``

        Parameters
        ----------
        group_name : str
            The name of the group
        sids : list of str
            The sample ids in the group

        Returns
        -------
        GroupResults
            The trajectory results for the given group

        Raises
        ------
        RuntimeError
            If sids is an empty list
        """
        # We multiply the coord values with the prop_expl
        trajectories = self._coords.loc[sids] * self._prop_expl

        if trajectories.empty:
            # Raising a RuntimeError since in a usual execution this should
            # never happen. The only way this can happen is if the user
            # directly calls this method, which shouldn't be done
            # (that's why the method is private)
            raise RuntimeError("No samples to process, an empty list cannot "
                               "be processed")

        # The weighting can only be done over trajectories with a length
        # greater than 1
        if self._weighted and len(sids) > 1:
            trajectories_copy = deepcopy(trajectories)
            try:
                trajectories = _weight_by_vector(trajectories_copy,
                                                 self._weighting_vector[sids])
            except (FloatingPointError, ValueError):
                self._message_buffer.append("Could not weight group, no "
                                            "gradient in the the "
                                            "weighting vector.\n")
                trajectories = trajectories_copy

        return self._compute_trajectories_results(group_name,
                                                  trajectories.loc[sids])

    def _compute_trajectories_results(self, group_name, trajectories):
        r"""Do the actual trajectories computation over trajectories

        Parameters
        ----------
        group_name : str
            The name of the group
        trajectories : pandas.DataFrame
            The sorted trajectories for each sample in the group

        Raises
        ------
        NotImplementedError
            This is the base class
        """
        raise NotImplementedError("No algorithm is implemented on the base "
                                  "class.")


class AverageGradientANOVA(GradientANOVA):
    r"""Perform trajectory analysis using the RMS average algorithm

    For each group in a category, it computes the average point among the
    samples in such group and then computes the norm of each sample from the
    averaged one.

    See Also
    --------
    GradientANOVA
    """

    _alg_name = 'avg'

    def _compute_trajectories_results(self, group_name, trajectories):
        r"""Do the actual trajectory computation over trajectories

        Parameters
        ----------
        group_name : str
            The name of the group
        trajectories : pandas.DataFrame
            The sorted trajectories for each sample in the group

        Returns
        -------
        GroupResults
            The trajectory results for `group_name` using the average
            trajectories method
        """
        center = np.average(trajectories, axis=0)
        if len(trajectories) == 1:
            trajectory = np.array([np.linalg.norm(center)])
            calc = {'avg': trajectory[0]}
        else:
            trajectory = np.array([np.linalg.norm(row[1].to_numpy() - center)
                                   for row in trajectories.iterrows()])
            calc = {'avg': np.average(trajectory)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return GroupResults(group_name, trajectory, np.mean(trajectory),
                            calc, msg)


class TrajectoryGradientANOVA(GradientANOVA):
    r"""Perform trajectory analysis using the RMS trajectory algorithm

    For each group in a category, each component of the result trajectory is
    computed as taking the sorted list of samples in the group and taking the
    norm of the coordinates of the 2nd sample minus 1st sample, 3rd sample
    minus 2nd sample and so on.

    See Also
    --------
    GradientANOVA
    """

    _alg_name = 'trajectory'

    def _compute_trajectories_results(self, group_name, trajectories):
        r"""Do the actual trajectory computation over trajectories

        Parameters
        ----------
        group_name : str
            The name of the group
        trajectories : pandas.DataFrame
            The sorted trajectories for each sample in the group

        Returns
        -------
        GroupResults
            The trajectory results for `group_name` using the trajectory
            method
        """
        if len(trajectories) == 1:
            trajectory = np.array([np.linalg.norm(trajectories)])
            calc = {'2-norm': trajectory[0]}
        else:
            # Loop through all the rows in trajectories and create '2-norm'
            # by taking the norm of the 2nd row - 1st row, 3rd row - 2nd row...
            trajectory = \
                np.array([np.linalg.norm(trajectories.iloc[i+1].to_numpy() -
                                         trajectories.iloc[i].to_numpy())
                          for i in range(len(trajectories) - 1)])
            calc = {'2-norm': np.linalg.norm(trajectory)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return GroupResults(group_name, trajectory, np.mean(trajectory),
                            calc, msg)


class FirstDifferenceGradientANOVA(GradientANOVA):
    r"""Perform trajectory analysis using the first difference algorithm

    It calculates the norm for all the time-points and then calculates the
    first difference for each resulting point

    See Also
    --------
    GradientANOVA
    """

    _alg_name = 'diff'

    def _compute_trajectories_results(self, group_name, trajectories):
        r"""Do the actual trajectory computation over trajectories

        Parameters
        ----------
        group_name : str
            The name of the group
        trajectories : pandas.DataFrame
            The sorted trajectories for each sample in the group

        Returns
        -------
        GroupResults
            The trajectory results for `group_name` using the first difference
            method
        """
        if len(trajectories) == 1:
            trajectory = np.array([np.linalg.norm(trajectories)])
            calc = {'mean': trajectory[0], 'std': 0}
        elif len(trajectories) == 2:
            trajectory = np.array([np.linalg.norm(trajectories[1] -
                                                  trajectories[0])])
            calc = {'mean': trajectory[0], 'std': 0}
        else:
            vec_norm = \
                np.array([np.linalg.norm(trajectories.iloc[i+1].to_numpy() -
                                         trajectories.iloc[i].to_numpy())
                          for i in range(len(trajectories) - 1)])
            trajectory = np.diff(vec_norm)
            calc = {'mean': np.mean(trajectory), 'std': np.std(trajectory)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return GroupResults(group_name, trajectory, np.mean(trajectory),
                            calc, msg)


class WindowDifferenceGradientANOVA(GradientANOVA):
    r"""Perform trajectory analysis using the modified first difference
    algorithm

    It calculates the norm for all the time-points and subtracts the mean of
    the next number of elements specified in `window_size` and the current
    element.

    Parameters
    ----------
    coords : pandas.DataFrame
        The coordinates for each sample id
    prop_expl : array like
        The numpy 1-D array with the proportion explained by each axis in
        coords
    metadata_map : pandas.DataFrame
        The metadata map, indexed by sample ids and columns are metadata
        categories
    window_size : int or long
        The window size to use while computing the differences

    Raises
    ------
    ValueError
        If the window_size is not a positive integer

    See Also
    --------
    GradientANOVA
    """

    _alg_name = 'wdiff'

    @experimental(as_of="0.4.0")
    def __init__(self, coords, prop_expl, metadata_map, window_size, **kwargs):
        super(WindowDifferenceGradientANOVA, self).__init__(coords, prop_expl,
                                                            metadata_map,
                                                            **kwargs)

        if not isinstance(window_size, Integral) or window_size < 1:
            raise ValueError("The window_size must be a positive integer")

        self._window_size = window_size

    def _compute_trajectories_results(self, group_name, trajectories):
        r"""Do the actual trajectory computation over trajectories

        If the first difference cannot be calculated of the provided window
        size, no difference is applied and a message is added to the results.

        Parameters
        ----------
        group_name : str
            The name of the group
        trajectories : pandas.DataFrame
            The sorted trajectories for each sample in the group

        Returns
        -------
        GroupResults
            The trajectory results for `group_name` using the windowed
            difference method
        """
        if len(trajectories) == 1:
            trajectory = np.array([np.linalg.norm(trajectories)])
            calc = {'mean': trajectory, 'std': 0}
        elif len(trajectories) == 2:
            trajectory = np.array([np.linalg.norm(trajectories[1] -
                                                  trajectories[0])])
            calc = {'mean': trajectory, 'std': 0}
        else:
            vec_norm = \
                np.array([np.linalg.norm(trajectories.iloc[i+1].to_numpy() -
                                         trajectories.iloc[i].to_numpy())
                          for i in range(len(trajectories) - 1)])
            # windowed first differences won't be able on every group,
            # specially given the variation of size that a trajectory tends
            # to have
            if len(vec_norm) <= self._window_size:
                trajectory = vec_norm
                self._message_buffer.append("Cannot calculate the first "
                                            "difference with a window of size "
                                            "(%d)." % self._window_size)
            else:
                # Replicate the last element as many times as required
                for idx in range(0, self._window_size):
                    vec_norm = np.append(vec_norm, vec_norm[-1:], axis=0)
                trajectory = []
                for idx in range(0, len(vec_norm) - self._window_size):
                    # Meas has to be over axis 0 so it handles arrays of arrays
                    element = np.mean(vec_norm[(idx + 1):
                                               (idx + 1 + self._window_size)],
                                      axis=0)
                    trajectory.append(element - vec_norm[idx])
                trajectory = np.array(trajectory)

            calc = {'mean': np.mean(trajectory), 'std': np.std(trajectory)}

        msg = ''.join(self._message_buffer) if self._message_buffer else None
        # Reset the message buffer
        self._message_buffer = []
        return GroupResults(group_name, trajectory, np.mean(trajectory),
                            calc, msg)
