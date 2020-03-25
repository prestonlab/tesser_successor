# creating and testing model RDMs
import scipy.spatial.distance as dist
import os
from glob import glob
import numpy as np

def rdm(matrix):
    """ Computes the representational dissimilarity matrix for one matrix. """
    return dist.squareform(dist.pdist(matrix, 'correlation'))


def multiple_rdm(SR_matrices):
    """ Computes the representational dissimilarity matrix for a list of matrices. """

    rdm_matrices = {}
    for part in [1, 2]:
        for run in range(1, 7):
            try:
                rdm_matrices[(part, run)] = rdm(SR_matrices[(part, run)])
            except KeyError:
                pass
    return rdm_matrices


def load_betas(data_dir, subject_num, roi): 
    """ Computes the representational dissimilarity matrix for one matrix. """
    
    roi_dir = os.path.join(data_dir, 'item_betas', 'roi')

    # look for directories with the correct pattern
    file_search = glob(os.path.join(roi_dir + '/' + roi, f'pattern_{subject_num}.txt'))
    if len(file_search) != 1:
        raise IOError(f'Problem finding ROI or subject directory for {roi} and {subject_num}')
    pattern_file = file_search[0]
    
    # read log, fixing problem with spaces in column names
    this_pattern = np.loadtxt(pattern_file)

    null_list = []
    for run in range(1, 7):
        this_run_len = run*151
        null_1 = this_run_len - 151
        null_list.append(null_1)
        null_2 = (this_run_len - 151) + 1
        null_list.append(null_2)
        null_3 = this_run_len - 2
        null_list.append(null_3)
        null_4 = this_run_len - 1
        null_list.append(null_4)
        
    # remove the fixation trials in the runs (these are just filler trials, i.e. the 4 null trials above)
    this_reformat_pattern = np.delete(this_pattern, null_list, axis=0)
    
    #return this_reformat_pattern
    return this_reformat_pattern
