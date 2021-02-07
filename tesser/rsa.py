"""Representational similarity analysis of objects after learning."""

from scipy.spatial import distance
import os
from glob import glob
import numpy as np
import pandas as pd
import scipy.spatial.distance as sd
from mindstorm import prsa
from tesser import util


def get_roi_sets():
    """Get a list of rois included in a set."""
    rois = {
        'hpc3': ['b_hip_ant', 'b_hip_body', 'b_hip_tail'],
        'mpfc9': ['10m', '10p', '10r', '11m', '14c', '14r', '24', '25', '32pl'],
    }
    return rois


def parse_rois(roi_list):
    """Parse an roi spec list."""
    roi_sets = get_roi_sets()
    split_rois = roi_list.split(',')
    rois = []
    for roi in split_rois:
        if roi in roi_sets:
            rois.extend(roi_sets[roi])
        else:
            rois.append(roi)
    return rois


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
    rsa_run_df = rsa_df[this[run - 1]: this[run], this[run - 1]: this[run]]
    return rsa_run_df


def load_zrep(data_dir, subject_num, roi):
    """ Computes the representational dissimilarity matrix for one matrix. """
    roi_dir = os.path.join(data_dir, 'item_zrep', 'roi')

    # look for directories with the correct pattern
    file_search = glob(os.path.join(roi_dir + '/' + roi, f'pattern_{subject_num}.txt'))
    if len(file_search) != 1:
        raise IOError(
            f'Problem finding ROI or subject directory for {roi} and {subject_num}'
        )
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
    return this_reformat_pattern


def load_vol_info(study_dir, subject):
    """Load volume information for all runs for a subject."""
    data_list = []
    columns = ['trial', 'onset', 'tr', 'sequence_type', 'trial_type', 'duration']
    runs = list(range(1, 7))
    for i, run in enumerate(runs):
        vol_file = os.path.join(
            study_dir, 'batch', 'analysis', 'rsa_beta', 'rsa_event_textfiles',
            f'tesser_{subject}_run{run}_info.txt'
        )
        run_data = pd.read_csv(vol_file, names=columns)
        run_data['duration'] = 1
        run_data['run'] = run
        data_list.append(run_data)
    data = pd.concat(data_list, axis=0)
    return data


def load_scram_info(event_dir):
    """Load event info for all runs."""

    # search for a file with the correct name formatting
    file_pattern = 'tesser_scram_volinfo.txt'
    file_search = glob(os.path.join(event_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding event info.')
    event_file = file_search[0]

    # read log, fixing problem with spaces in column names
    event_info = pd.read_csv(event_file, sep='\,', skipinitialspace=True, header=None)

    # name the columns:
    event_info = event_info.rename({0: 'run', 1: 'item', 2: 'comm', 3: 'bound'}, axis=1)
    return event_info


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
    """
    Pairs where conditions are equal.

    e.g. Pairs are from the same community for all three communities.
    """
    return x[:, None] == x[:, None].T


def pair_neq(x):
    """
    Pairs where conditions are not equal.

    e.g. Pairs are not selected from the same run.
    """
    return x[:, None] != x[:, None].T


def pair_and(x):
    """Pairs where conditions are both true."""
    return x[:, None] & x[:, None].T


# getting averages of matrix data, making sure that it is below diagonal etc.
def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""
    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat


def load_roi_brsa(res_dir, rois, subjects=None):
    """Load correlation matrices from BRSA results."""
    if subjects is None:
        subjects = util.subj_list()
    rdms = {
        roi: [
            np.load(os.path.join(res_dir, roi, f'sub-{subject}_brsa.npz'))['C']
            for subject in subjects
        ] for roi in rois
    }
    return rdms


def load_roi_mean_brsa(res_dir, rois, subjects=None):
    if subjects is None:
        subjects = util.subj_list()
    rdms = load_roi_brsa(res_dir, rois, subjects)
    n_subj = len(subjects)
    mrdm = {
        roi: np.mean(
            np.dstack([rdms[roi][subject] for subject in range(n_subj)]), axis=2
        ) for roi in rois
    }
    return mrdm


def load_roi_prsa(res_dir, roi, subjects=None):
    """Load z-statistic from permutation test results."""
    if subjects is None:
        subjects = util.subj_list()

    z = []
    for subject in subjects:
        subj_file = os.path.join(res_dir, roi, f'zstat_{subject}.csv')
        zdf = pd.read_csv(subj_file, index_col=0).T
        zdf.index = [subject]
        z.append(zdf)
    df = pd.concat(z)
    return df


def load_net_prsa(rsa_dir, block, model_set, rois, subjects=None):
    """Load z-statistics for a model for a set of ROIs."""
    # get the directory with results for this model set and block
    res_dir = os.path.join(rsa_dir, f'prsa_{block}_{model_set}')
    if not os.path.exists(res_dir):
        raise IOError(f'Results directory not found: {res_dir}')

    # load each ROI
    df_list = []
    for roi in rois:
        rdf = load_roi_prsa(res_dir, roi, subjects)
        mdf = pd.DataFrame({'subject': rdf.index, 'roi': roi}, index=rdf.index)
        full = pd.concat((mdf, rdf), axis=1)

        df_list.append(full)
    df = pd.concat(df_list, ignore_index=True)
    return df


def net_prsa_perm(df, model, n_perm=1000, beta=0.05):
    """Test ROI correlations using a permutation test."""
    # shape into matrix format
    rois = df.roi.unique()
    mat = df.pivot(index='subject', columns='roi', values=model)
    mat = mat.reindex(columns=rois)

    # run sign flipping test
    results = prsa.sign_perm(mat.to_numpy(), n_perm, beta)
    results.index = mat.columns
    return results
