"""Representational similarity analysis of objects after learning."""

from scipy.spatial import distance
import os
from glob import glob
import numpy as np
import scipy.spatial.distance as sd


def rdm(matrix):
    """ Computes the representational dissimilarity matrix for one matrix. """
    return distance.squareform(distance.pdist(matrix, 'correlation'))


def multiple_rdm(SR_matrices):
    """Representational dissimilarity matrices for a list of SR matrices."""

    rdm_matrices = {}
    for part in [1, 2]:
        for run in range(1, 7):
            try:
                rdm_matrices[(part, run)] = rdm(SR_matrices[(part, run)])
            except KeyError:
                pass
    return rdm_matrices


def load_rsa(data_dir, subject, roi):
    """Load RSA data by subject."""

    # search for a file with the correct name formatting
    file_pattern = f'{subject}_betas_{roi}.txt'
    file_search = glob(os.path.join(data_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject}.')
    rsa_file = file_search[0]

    # read log, fixing problem with spaces in column names
    rsa_df = np.loadtxt(rsa_file)

    return rsa_df


def rsa_run(rsa_df, run):
    """Load RSA data by run."""
    this = list(range(0, 1057, 151))
    rsa_run_df = rsa_df[this[run - 1]:this[run], this[run - 1]:this[run]]

    return rsa_run_df


def load_betas(data_dir, subject_num, roi):
    """ Computes the representational dissimilarity matrix for one matrix. """

    roi_dir = os.path.join(data_dir, 'item_betas', 'roi')

    # look for directories with the correct pattern
    file_search = glob(
        os.path.join(roi_dir + '/' + roi, f'pattern_{subject_num}.txt'))
    if len(file_search) != 1:
        raise IOError(f'Problem finding ROI or subject directory for '
                      f'{roi} and {subject_num}')
    pattern_file = file_search[0]

    # read log, fixing problem with spaces in column names
    this_pattern = np.loadtxt(pattern_file)

    null_list = []
    for run in range(1, 7):
        this_run_len = run * 151
        null_1 = this_run_len - 151
        null_list.append(null_1)
        null_2 = (this_run_len - 151) + 1
        null_list.append(null_2)
        null_3 = this_run_len - 2
        null_list.append(null_3)
        null_4 = this_run_len - 1
        null_list.append(null_4)

    # remove the fixation trials in the runs
    # (these are just filler trials, i.e. the 4 null trials above)
    this_reformat_pattern = np.delete(this_pattern, null_list, axis=0)

    #  return this_reformat_pattern
    return this_reformat_pattern


def exclude_rsa(rsa, exclude_n):
    for d in range(0, exclude_n):
        diag_pos = np.diagonal(rsa, d)
        diag_neg = np.diagonal(rsa, -d)
        # It's not writable. MAKE it writable.
        diag_pos.setflags(write=True)
        diag_neg.setflags(write=True)
        diag_pos.fill('NaN')
        diag_neg.fill('NaN')
    return rsa


def pair_eq(x):
    """Pairs where conditions are equal."""
    """e.g. Pairs are from the same community for all three communities."""

    return x[:, None] == x[:, None].T


def pair_neq(x):
    """Pairs where conditions are not equal."""
    """e.g. Pairs are not selected from the same run."""
    return x[:, None] != x[:, None].T


def pair_and(x):
    """Pairs where conditions are both true."""
    """e.g. Pairs ...??"""
    return x[:, None] & x[:, None].T


# getting averages of matrix data, making sure that it is below diagonal etc.
def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""

    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat
