# -*- coding: utf-8 -*-
'''
Implementation of the Knox test for space-time clustering in point data, using a
permutation-based test as used in Johnson et al. (2007) 'Space–Time Patterns of 
Risk: A Cross National Assessment of Residential Burglary Victimization'.

Author: Toby Davies
'''


import numpy as np 
import pandas as pd
import scipy.stats as ss
from sklearn.neighbors import radius_neighbors_graph


def compute_pairs(xy, t, s_bands, t_bands, metric="euclidean", interval_side="left"):
    # Get inter-event distances up to maximum required threshold
    if isinstance(xy, pd.DataFrame):
        s_pairs = radius_neighbors_graph(xy.values, s_bands[-1], mode="distance", metric=metric)
    else:
        s_pairs = radius_neighbors_graph(xy, s_bands[-1], mode="distance", metric=metric)
    if isinstance(t, pd.Series):
        t_pairs = radius_neighbors_graph(t.values.reshape(-1,1), t_bands[-1], mode="distance")
    else:
        t_pairs = radius_neighbors_graph(t.reshape(-1,1), t_bands[-1], mode="distance")
    # Convert distances to bands
    # Add 1 to make them non-zero
    s_pairs.data = np.searchsorted(s_bands, s_pairs.data, interval_side) + 1
    t_pairs.data = np.searchsorted(t_bands, t_pairs.data, interval_side) + 1
    
    if interval_side=="right":
        t_pairs.data[t_pairs.data == (t_bands.shape[0] + 1)] = 0
        t_pairs.eliminate_zeros()
    
    return s_pairs, t_pairs


def calculate_significance(conting_observed, conting_perm):
    
    #Stack the observed contingency below the permutation ones
    conting_all = np.dstack((conting_observed, conting_perm))
    
    #Get the sorted indices depth-wise and find the position of the observed
    indices = np.argsort(conting_all, axis=-1)
    ranks = np.argmax(indices == 0, axis=-1)

    #Calculate p-value by referring to the ranks of the observed values
    p_values = 1 - (ranks / np.float64(conting_all.shape[-1]))
    
    #For each contingency cell, get the median value across all permutation iterations
    medians_perm = np.median(conting_perm, axis=-1)
    means_perm = np.mean(conting_perm, axis=-1)
    
    #Knox ratios found by dividing observed counts by medians
    ratios_median = conting_observed / medians_perm
    ratios_mean = conting_observed / means_perm
    
    #Also calculate z-scores of observed values as alternative measure of deviation
    z_scores = ss.zscore(conting_all, axis=-1)[...,0]
    
    return p_values, ratios_median, ratios_mean, z_scores


def compute_contingency(s_pairs, t_pairs, s_bands, t_bands):
    st_pairs = s_pairs.multiply(t_pairs)
    st_pairs_r, st_pairs_c = st_pairs.nonzero()
    
    s_bins = s_pairs[st_pairs_r, st_pairs_c].A1 - 1
    t_bins = t_pairs[st_pairs_r, st_pairs_c].A1 - 1
    
    bin_edges = [np.arange(s_bands.shape[0]+1), np.arange(t_bands.shape[0]+1)]
    conting, _xedges, _yedges  = np.histogram2d(s_bins, t_bins, bin_edges)
    
    return conting / 2.0


def knox_test(xy, t, s_bands, t_bands, n_iter, metric="euclidean", interval_side="left",
              seed=None):
    """
    Perform the Knox test
    
    Parameters
    ----------
    
    :param xy: Array of shape (N, 2) representing the x and y coordinates of the 
    event location in metric units
    
    :param t: 1-dimensional array of length n_events representing timestamps of events,
    expressed in numerical form representing some unit (e.g. days)
    
    :param s_bands: Upper limits of spatial thresholds
    
    :param t_bands: Upper limits of temporal thresholds
    
    :param n_iter: Number of iterations to perform
    
    :param metric: Distance metric for spatial computations - currently supported options
    "euclidean" and "manhattan"
    
    :param interval_side: Determines the side on which the band intervals are open -
    default is "left"
    
    :param seed: random seed to initialise test
    """
    if interval_side not in ('left', 'right'):
        raise ValueError("Supported values for interval_side are 'left', 'right'.")
    if xy.shape[0]!=t.shape[0]:
        raise ValueError("xy and t must have the same number of points.")
    # Set the random seed if required
    if seed is not None:
        np.random.seed(seed)
    # Convert bands to arrays
    s_bands = np.asarray(s_bands)
    t_bands = np.asarray(t_bands)
    # Calculate thresholded pairwise distances in each dimension
    s_pairs, t_pairs = compute_pairs(xy, t, s_bands, t_bands, metric, interval_side)
    # Calculate observed contingency
    conting_observed = compute_contingency(s_pairs, t_pairs, s_bands, t_bands)
    
    # Initialise contingency array for permutation iteraions
    s_shape = s_bands.shape[0]
    t_shape = t_bands.shape[0]
    conting_perm = np.zeros((s_shape, t_shape, n_iter), dtype=np.float64)
    
    # Define index for shuffling
    perm = np.arange(xy.shape[0])
    for k in range(n_iter):
        # Shuffle the indices
        np.random.shuffle(perm)
        perm_inv = np.argsort(perm)
        # Apply shuffle to spatial distance matrix (in CSR form)
        s_pairs_shuff = s_pairs[perm, :]
        s_pairs_shuff.indices = perm_inv[s_pairs_shuff.indices]
        # Compute contingency for shuffled data and add to master array
        this_conting_perm = compute_contingency(s_pairs_shuff, t_pairs, s_bands, t_bands)
        conting_perm[:,:,k] = this_conting_perm
    
    # Calculate statistical outputs
    p_values, ratios_median, ratios_mean, z_scores = calculate_significance(conting_observed, conting_perm)
    
    result = KnoxResult(conting_observed, p_values, ratios_median, ratios_mean,
                        z_scores, s_bands, t_bands, n_iter, metric, interval_side)
    
    return result


class KnoxResult(object):
    
    def __init__(self, conting_observed, p_values, ratios_median, ratios_mean,
                 z_scores, s_bands, t_bands, n_iter, metric, interval_side):
        
        self.s_bands = s_bands
        self.t_bands = t_bands
        self.interval_side = interval_side
        
        if self.interval_side == 'left':
            left_bracket, right_bracket = "(", "]"
        else:
            left_bracket, right_bracket = "[", ")"
        
        self.s_labels = ["{}{}, {}{}".format(left_bracket, l, r, right_bracket) for l, r in zip(self.s_bands, self.s_bands[1:])]
        self.s_labels.insert(0, "[0, {}{}".format(self.s_bands[0], right_bracket))
        self.t_labels = ["{}{}, {}{}".format(left_bracket, l, r, right_bracket) for l, r in zip(self.t_bands, self.t_bands[1:])]
        self.t_labels.insert(0, "[0, {}{}".format(self.t_bands[0], right_bracket))
        
        self.conting_observed = pd.DataFrame(conting_observed, index=self.s_labels, columns=self.t_labels)
        self.p_values = pd.DataFrame(p_values, index=self.s_labels, columns=self.t_labels)
        self.ratios_median = pd.DataFrame(ratios_median, index=self.s_labels, columns=self.t_labels)
        self.ratios_mean = pd.DataFrame(ratios_mean, index=self.s_labels, columns=self.t_labels)
        self.z_scores = pd.DataFrame(z_scores, index=self.s_labels, columns=self.t_labels)
        self.n_iter = n_iter
        self.metric = metric
        
    
    def plot_ratios_median(self, vmin=0, vmax=2, cmap="autumn_r"):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        cm = plt.get_cmap(cmap)
        norm = colors.Normalize(vmin, vmax)
        mapper = lambda x: 'background-color: {}'.format(colors.to_hex(cm(norm(x))))
        return self.ratios_median.style.applymap(mapper)
        
        