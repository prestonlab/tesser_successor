# creating and testing model RDMs
import scipy as sp
import scipy.spatial.distance as dist


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

