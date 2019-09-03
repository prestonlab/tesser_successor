# creating and testing model RDMs
import scipy as sp
import scipy.spatial.distance as dist


def rdm(matrix):
    """ Computes the representational dissimilarity matrix for one matrix. """
    return dist.squareform(dist.pdist(matrix, 'correlation'))


def multiple_rdm(matrices):
    """ Computes the representational dissimilarity matrix for a list of matrices. """
    rdm_matrices = []
    for matrix in matrices:
        rdm_matrices.append(rdm(matrix))
    return rdm_matrices
