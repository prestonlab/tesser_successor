import pytest
import numpy as np
import pandas as pd

@pytest.fixture
def matrix_counts():
    expected = np.array(
        [[11., 8., 11.],
         [11., 11., 9.],
         [8., 12., 17.]]
    )
    return expected

def test_transition_counts(matrix_counts):
    objnum = np.array([1, 0, 1, 0, 1, 0, 0, 2, 0, 2, 1, 0, 2, 2, 1, 2, 0, 0, 0, 0, 0, 2,
       2, 2, 0, 2, 2, 0, 2, 2, 2, 1, 2, 1, 1, 0, 2, 2, 1, 2, 2, 2, 1, 1,
       1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 2, 2, 2, 0, 0, 0, 2, 2, 2, 0, 0, 1,
       0, 1, 2, 0, 0, 1, 1, 0, 1, 2, 2, 1, 0, 2, 1, 1, 1, 0, 0, 1, 1, 2,
       2, 2, 1, 0, 0, 2, 2, 1, 1, 1, 1])
    matrix = np.zeros([3, 3])
    for i, j in enumerate(objnum):
        try:
            matrix[j,objnum[i+1]]+=1 # did not have to subtract minus 1
        except:
            KeyError
    np.testing.assert_allclose(matrix, matrix_counts)
    

@pytest.fixture
def final_matrix():
    expected = np.array(
      [[0.36666667, 0.26666667, 0.36666667],
       [0.35483871, 0.35483871, 0.29032258],
       [0.21621622, 0.32432432, 0.45945946]]
    )
    return expected
    
def test_transition_matrix(final_matrix):
    objnum = np.array([1, 0, 1, 0, 1, 0, 0, 2, 0, 2, 1, 0, 2, 2, 1, 2, 0, 0, 0, 0, 0, 2,
       2, 2, 0, 2, 2, 0, 2, 2, 2, 1, 2, 1, 1, 0, 2, 2, 1, 2, 2, 2, 1, 1,
       1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 2, 2, 2, 0, 0, 0, 2, 2, 2, 0, 0, 1,
       0, 1, 2, 0, 0, 1, 1, 0, 1, 2, 2, 1, 0, 2, 1, 1, 1, 0, 0, 1, 1, 2,
       2, 2, 1, 0, 0, 2, 2, 1, 1, 1, 1])
    matrix = np.zeros([3, 3])
    for i, j in enumerate(objnum):
        try:
            matrix[j,objnum[i+1]]+=1  # did not have to subtract minus 1
        except:
            KeyError
    for row in range(3):
        matrix[row] /= np.sum(matrix[row])
    np.testing.assert_allclose(matrix, final_matrix)
    
@pytest.fixture
def matrix_counts_with_missing_obj():
    expected = np.array(
      [[ 0.,  0.,  0.],
       [ 0., 27., 25.],
       [ 0., 26., 20.]]
    )
    return expected

def test_transition_counts(matrix_counts_with_missing_obj):
    objnum = np.array(
       [2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2,
       1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 2, 2,
       1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 1,
       2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1,
       2, 1, 2, 2, 1, 1, 1, 2, 1, 2, 1]
    )
    
    matrix = np.zeros([3, 3])
    for i, j in enumerate(objnum):
        try:
            matrix[j,objnum[i+1]]+=1  # did not have to subtract minus 1
        except:
            KeyError
    np.testing.assert_allclose(matrix, matrix_counts_with_missing_obj)
    

    
@pytest.fixture
def final_matrix_counts_with_missing_obj():
    expected = np.array(
      [[np.nan    , np.nan    ,    np. nan],
       [0.        , 0.51923077, 0.48076923],
       [0.        , 0.56521739, 0.43478261]]
    )
    return expected

def test_transition_error(final_matrix_counts_with_missing_obj):
    objnum = np.array(
       [2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2,
       1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1, 2, 2, 1, 2, 2,
       1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 1,
       2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 1, 1, 1,
       2, 1, 2, 2, 1, 1, 1, 2, 1, 2, 1]
    )

    matrix = np.zeros([3, 3])
    for i, j in enumerate(objnum):
        try:
            matrix[j,objnum[i+1]]+=1  # did not have to subtract minus 1
        except:
            KeyError
    for row in range(3):
        matrix[row] = np.divide(matrix[row],np.sum(matrix[row]))
    np.testing.assert_allclose(matrix, final_matrix_counts_with_missing_obj)