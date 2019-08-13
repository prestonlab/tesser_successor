import scipy as sp
import scipy.spatial.distance as dist

def rda (matrix):
    ''' Computes the representational dissimilarity matrix for one matrix. '''
    return dist.squareform (dist.pdist (matrix, 'correlation'))

def multiple_rda (matrices):
    ''' Computes the representational dissimilarity matrix for a list of matrices. '''
    rda_matrices = []
    for matrix in matrices:
        rda_matrices.append (rda (matrix))
    return rda_matrices