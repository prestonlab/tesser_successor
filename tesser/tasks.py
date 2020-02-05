"""Module for working with grouping data."""

import numpy as np
import scipy.spatial.distance as sd


#make a function that gets the average inference task data for a participant

#make a function that gets the average inference task data for all participants and
#then does a median split

def group_dist(group_mat, distance='Euclidean'):
    """Calculate object distance from grouping task."""

    rows, cols = np.nonzero(group_mat)
    objects = group_mat[rows, cols]
    coord = np.vstack((rows, cols)).T
    sort_ind = np.argsort(objects)
    coord_sort = coord[sort_ind]
    group_dist_mat = sd.squareform(sd.pdist(coord_sort, distance))
    return group_dist_mat


def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""

    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat
