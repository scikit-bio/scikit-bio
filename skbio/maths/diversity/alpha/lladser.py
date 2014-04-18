#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from .base import _validate


# NOT TESTED: NEED TEST DATA!
def lladser_pe(counts, r=10, **args):
    """Single point estimate of the conditional uncovered probability

    This function is just a wrapper around the full point estimator,
    intended to be called fo a single best estimate on a complete sample.
    """
    counts = _validate(counts)
    sample = _expand_counts(counts)
    np.random.shuffle(sample)
    try:
        pe = list(_lladser_point_estimates(sample, r))[-1][0]
    except IndexError:
        pe = 'NaN'
    return pe


# NOT TESTED: NEED TEST DATA!
def lladser_ci(counts, r, alpha=0.95, f=10, ci_type='ULCL', **args):
    """Single CI of the conditional uncovered probability

    This function is just a wrapper around the full point estimator,
    intended to be called for a single best estimate on a complete sample.
    """
    counts = _validate(counts)
    sample = _expand_counts(counts)
    np.random.shuffle(sample)
    try:
        pe = list(_lladser_ci_series(sample, r))[-1]
    except IndexError:
        pe = ('NaN', 'NaN')
    return pe


def _expand_counts(counts):
    """Converts vector of counts at each index to vector of indices."""
    result = []
    for i, c in enumerate(counts):
        result.append(np.zeros(c, int) + i)
    return np.concatenate(result)


def _lladser_point_estimates(sample, r=10):
    """Series of point estimates of the conditional uncovered probability

    sample: series of random observations
    r: Number of new colors that are required for the next prediction

    This is the point estimator described in Theorem 2 (i):
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.
    Returns: Each new color yields 3-tuple:
         - point estimate
         - position in sample of prediction
         - random variable from poisson process (mostly to make testing easier)
    """
    if(r <= 2):
        raise ValueError("r must be >=3")
    for count, seen, cost, i in _get_interval_for_r_new_species(sample, r):
        t = np.random.gamma(count, 1)
        point_est = (r - 1) / t
        yield(point_est, i, t)


def _get_interval_for_r_new_species(seq, r):
    """For seq of samples compute interval between r new species.

    seq: series of observations (the actual sample, not the frequencies)
    r: Number of new colors to that need to be observed for a new interval

    Imagine an urn with colored balls. Given a drawing of balls from the urn,
    compute how many balls need to be looked at to discover r new colors.
    Colors can be repeated.

    Returns: for each new color seen for the first time, yield:
             - length of interval, i.e. number of observations looked at
             - the set of seen colors
             - position in seq after seeing the last new color (end of interval)
             - position in seq where interval is started
    """
    seen = set()
    seq_len = len(seq)
    for i in range(seq_len):
        curr = seq[i]  # note: first iteration is after looking at first char
        # bail out if there's nothing new
        if curr in seen:
            continue
        else:
            seen.add(curr)

        # otherwise, need to see distance to get k colors
        unseen = 0
        j = i + 1
        while unseen < r and j < seq_len:
            if seq[j] not in seen:
                unseen += 1
            j += 1  # note: increments after termination condition

        count = j - i - 1  # the interval to see r new colors
        cost = j  # the position in seq after seeing r new ones
        if (not count) or (unseen < r):  # bail out if not enough unseen
            raise StopIteration
        yield count, set(seen), cost, i


def _lladser_ci_series(seq, r, alpha=0.95, f=10, ci_type='ULCL'):
    """Construct r-color confidence intervals for uncovered conditional prob.

    seq: Input is a sequence of colors (the actual sample, not the counts)
    r  : Number of new colors that are required for the next prediction
    alpha: desired confidence level
    f: ratio between upper and lower bound
    ci_type: type of confidence interval. One of:
             ULCL: upper and lower bounds with conservative lower bound
             ULCU: upper and lower woth conservative upper bound
             U: Upper bound only, lower bound fixed to 0
             L: Lower bound only, upper bound fixed to 1

    Returns: One CI prediction for each new color that is detected and where.
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via
    Poissonization: Accurate Measurements of the Microbial Unknown" PLoS 2011.

    """
    for count, seen, cost, i in _get_interval_for_r_new_species(seq, r):
        t = np.random.gamma(count, 1)
        yield _lladser_ci_from_r(r, t, alpha, f, ci_type)


def _lladser_ci_from_r(r, t, alpha=0.95, f=10, ci_type='ULCL'):
    """Construct r-color confidence intervals for uncovered conditional prob.

    r: Number of new colors that are required for the next prediction
    t: A value from the gamma distribution gamma (count,1)
    alpha: desired confidence level
    f: ratio between upper and lower bound
    ci_type: type of confidence interval. One of:
             ULCL: upper and lower bounds with conservative lower bound
             ULCU: upper and lower woth conservative upper bound
             U: Upper bound only, lower bound fixed to 0
             L: Lower bound only, upper bound fixed to 1

    This is the formula that is described in Theorem 2 iii
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.

    Returns: A confidence interval that contains the true conditional
             uncovered probability with a probability of 100% * alpha.
    """
    if ci_type == 'U':
        upper_bound = _upper_confidence_bound(r, alpha) / t
        return (0, upper_bound)
    elif ci_type == 'L':
        lower_bound = _lower_confidence_bound(r, alpha) / t
        return (lower_bound, 1)

    bound_params = _ul_confidence_bounds(f, r, alpha)
    if ci_type == 'ULCL':
        bound_param = bound_params[0]
    elif ci_type == 'ULCU':
        bound_param = bound_params[1]
    else:
        raise ValueError("Unknown ci_type: %s" % (ci_type))

    upper_bound = bound_param * f / t
    lower_bound = bound_param / t

    # make sure upper bound is at most 1
    if (upper_bound > 1):
        upper_bound = 1

    return lower_bound, upper_bound


def _upper_confidence_bound(r, alpha):
    """Get constant for confidence interval with lower bound fixed to 0

    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    Compute constant b according to Theorem 2 iii with a=0
    aka c_0 from Table 3
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via
    Poissonization: Accurate Measurements of the Microbial Unknown" PLoS 2011.

    Returns: Constant c such that the confidence interval is [0,c/T_r]

    """
    alpha = round(alpha, 2)
    if not (alpha == 0.95):
        raise ValueError("alpha must be 0.95")

    data = {1: 2.995732274,
            2: 4.743864518,
            3: 6.295793622,
            4: 7.753656528,
            5: 9.153519027,
            6: 10.51303491,
            7: 11.84239565,
            8: 13.14811380,
            9: 14.43464972,
            10: 15.70521642,
            11: 16.96221924,
            12: 18.20751425,
            13: 19.44256933,
            14: 20.66856908,
            15: 21.88648591,
            16: 23.09712976,
            17: 24.30118368,
            18: 25.49923008,
            19: 26.69177031,
            20: 27.87923964,
            21: 29.06201884,
            22: 30.24044329,
            23: 31.41481021,
            24: 32.58538445,
            25: 33.75240327,
            50: 62.17105670,
            }
    try:
        return (data[r])
    except KeyError:
        raise ValueError("r must be between 1,..,25 or 50")


def _lower_confidence_bound(r, alpha):
    """Get constant for confidence interval with upper bound fixed to 1

    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    Compute constant b according to Theorem 2 iii with b=1
    aka c_3 from Table 3
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via
    Poissonization: Accurate Measurements of the Microbial Unknown" PLoS 2011.
    Returns: Constant c such that the confidence interval is [c/T_r, 1]

    """
    alpha = round(alpha, 2)
    if not (alpha == 0.95):
        raise ValueError("alpha must be 0.95")
    data = {1: 0.051293294,
            2: 0.355361510,
            3: 0.817691447,
            4: 1.366318397,
            5: 1.970149568,
            6: 2.613014744,
            7: 3.285315692,
            8: 3.980822786,
            9: 4.695227540,
            10: 5.425405697,
            11: 6.169007289,
            12: 6.924212514,
            13: 7.689578292,
            14: 8.463937522,
            15: 9.246330491,
            16: 10.03595673,
            17: 10.83214036,
            18: 11.63430451,
            19: 12.44195219,
            20: 13.25465160,
            21: 14.07202475,
            22: 14.89373854,
            23: 15.71949763,
            24: 16.54903871,
            25: 17.38212584
            }
    try:
        return (data[r])
    except KeyError:
        raise ValueError("r must be between 1,..,25 or 50")


def _ul_confidence_bounds(f, r, alpha):
    """returns confidence bounds based on ratio f and alpha

    f: desired ratio of upper to lower bound
    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    This function is just a lookup of some precomputed values.

    Returns: Constant c_1,c_2 such that the confidence interval is:
             [c_1/T_r, c_1*f/T_r] for conservative lower bound intervals and
             [c_2/T_r, c_2*f/T_r] for conservative upper bound intervals
    """
    alpha = round(alpha, 2)
    a = None
    b = None

    if (f, r, alpha) in PRECOMPUTED_TABLE:
        return PRECOMPUTED_TABLE[(f, r, alpha)]

    # all others combination are only computed for f=10
    # and alpha = 0.90, 0.95 and 0.99
    if f == 10 and r <= 50:
        if alpha in CBS and r < len(CBS[alpha]):
            (a, b) = CBS[alpha][r]

    if a is None or b is None:
        raise ValueError("No constants are precomputed for the combination of "
                         "f:%f, r:%d, and alpha:%.2f" % (f, r, alpha))

    return (a, b)


# Hack in some special values we used for the paper.
# Since Manuel needs to compute those semi-automatically
# using Maple, we pre-calculate only a few common ones

# precomputed table is {(f, r, alpha):(c_1, c_2)}
PRECOMPUTED_TABLE = {
    (2, 50, 0.95): (31.13026306, 38.94718565),
    (2, 33, 0.95): (22.3203508, 23.4487304),
    (1.5, 100, 0.95): (79.0424349, 83.22790086),
    (1.5, 94, 0.95): (75.9077267, 76.5492088),
    (2.5, 19, 0.95): (11.26109001, 11.96814857),

    # In the next block for each f, we report the smallest possible value
    # of r from table 4 in the paper
    (80, 2, 0.95): (0.0598276655, 0.355361510),
    (48, 2, 0.95): (0.1013728884, 0.355358676),
    (40, 2, 0.95): (0.1231379857, 0.355320458),
    (24, 2, 0.95): (0.226833483, 0.346045204),
    (20, 3, 0.95): (0.320984257, 0.817610455),
    (12, 3, 0.95): (0.590243030, 0.787721610),
    (10, 4, 0.95): (0.806026244, 1.360288674),
    (6, 6, 0.95): (1.8207383, 2.58658608),
    (5, 7, 0.95): (2.48303930, 3.22806682),
    (3, 14, 0.95): (7.17185045, 8.27008349),
    (2.5, 19, 0.95): (11.26109001, 11.96814857),
    (1.5, 94, 0.95): (75.9077267, 76.5492088),
    (1.25, 309, 0.95): (275.661191, 275.949782)
}


######
# Below are the values used for Theorem 3 iii
# Values hand computed by Manuel Lladser using Maple

CB_90 = [
    (None, None),  # 0, makes indexing easier
    (None, None),  # no feasible solution
    (None, None),  # no feasible solution
    (.5635941995, 1.095834700),
    (.6764656264, 1.744588615),
    (.8018565594, 2.432587343),
    (.9282215025, 3.151897973),
    (1.053433716, 3.894766804),
    (1.177158858, 4.656118177),
    (1.299491033, 5.432468058),
    (1.420604842, 6.221304605),  # 10
    (1.540665805, 7.020746595),
    (1.659812701, 7.829342026),
    (1.778158703, 8.645942495),
    (1.895796167, 9.469621185),
    (2.012801198, 10.29961731),
    (2.129237257, 11.13529724),
    (2.245157877, 11.97612664),
    (2.360608695, 12.82164994),
    (2.475628991, 13.67147502),
    (2.590252861, 14.52526147),  # 20
    (2.704510123, 15.38271151),
    (2.818427036, 16.24356290),
    (2.932026869, 17.10758326),
    (3.045330351, 17.97456551),
    (3.158356050, 18.84432420),
    (None, None),  # not computed
    (None, None),
    (None, None),
    (None, None),
    (3.719850286, 23.22944415),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.828910181, 32.13892224),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.924900191, 41.17906791)  # 50
]

CB_95 = [
    (None, None),  # 0
    (None, None),
    (None, None),
    (None, None),
    (.8060262438, 1.360288674),  # 4
    (.9240311584, 1.969902537),
    (1.053998892, 2.613007253),
    (1.185086998, 3.285315518),
    (1.315076337, 3.980822783),
    (4.695227540, 4.695227541),
    (1.570546801, 5.425405698),  # 10
    (1.696229569, 6.169007289),
    (1.820753729, 6.924212513),
    (1.944257622, 7.689578291),
    (2.066857113, 8.463937522),
    (2.188648652, 9.246330491),
    (2.309712994, 10.03595673),
    (2.430118373, 10.83214036),
    (2.549923010, 11.63430451),
    (2.669177032, 12.44195219),
    (2.787923964, 13.25465160),  # 20
    (2.906201884, 14.07202475),
    (3.024044329, 14.89373854),
    (3.141481020, 15.71949763),
    (3.258538445, 16.54903871),
    (3.375240327, 17.38212584),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (3.954097220, 21.59397923),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.093973695, 30.19573919),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.217105673, 38.96473258)  # 50
]

CB_99 = [
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (1.360316290, 1.768978323),
    (1.470856924, 2.329171347),
    (1.604478487, 2.906049304),
    (1.741759456, 3.507452949),
    (1.878809285, 4.130199076),  # 10
    (2.014632329, 4.771246173),
    (2.149044735, 5.428180734),
    (2.282101533, 6.099073460),
    (2.413917374, 6.782354878),
    (2.544610844, 7.476728267),
    (2.674289153, 8.181107778),
    (2.803045614, 8.894573463),
    (2.930960779, 9.616337916),
    (3.058104355, 10.34572103),
    (3.184536992, 11.08213063),  # 20
    (3.310311816, 11.82504734),
    (3.435475649, 12.57401269),
    (3.560070013, 13.32861956),
    (3.684131925, 14.08850448),
    (3.807694563, 14.85334135),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.41897094, 18.7424258),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.61643962, 26.7700386),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.79033616, 35.0324474)  # 50
]

CBS = {0.90: CB_90,
       0.95: CB_95,
       0.99: CB_99}
