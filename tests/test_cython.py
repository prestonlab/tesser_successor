import pytest
import cython
import numpy as np
import pandas as pd
from tesser import csr
from tesser import cfit
from tesser import util

@pytest.fixture
def cy_SR():
    expected =np.array(
        [[0.140625 , 0.1640625, 0.78125  ],
         [0.21875  , 0.4375   , 0.5      ],
         [0.5625   , 0.40625  , 0.125    ]]
    )

    return expected
    
def test_cy_learn(cy_SR):
    M = [[0,0,0],
         [0,0,0],
         [0,0,0]]
    M = np.double(M)
    onehot= [[1,0,0],
             [0,1,0],
             [0,0,1]]
    onehot = np.int32(onehot)
    envstep = np.array([2,1,0,2,1,1,2,0,2])
    envstep = envstep.astype(np.dtype('i'))
    alpha = 0.5
    gamma = 0.5
    n_states = 3
    csr.SR(envstep, gamma, alpha, M, n_states, onehot)
    
    np.testing.assert_allclose(M, cy_SR)
    
def dummy_induct_df():
    # construct a data frame with dummy data
    data ={}
    data['subjnum'] =np.hstack((np.zeros([1,33]),np.ones([1,33]),np.ones([1,33])+1))[0].astype('int') +100
    data['cue'] =np.random.randint(0,3,99)
    data['opt1'] = np.random.randint(0,3,99)
    data['opt2'] =np.random.randint(0,3,99)
    df = pd.DataFrame(data)

    return df

@pytest.fixture
def transition_matrix():
    expected = np.array(
      [[np.nan    , np.nan    ,    np. nan],
       [0.        , 0.51923077, 0.48076923],
       [0.        , 0.56521739, 0.43478261]]
    )
    return expected

def test_prob_choice(transition_matrix): # test what happens for nans
    sr = np.array(
        [[0.140625 , 0.1640625, 0.  ],
         [0.       , 0.       , 0.  ],
         [0.5625   , 0.40625  , 0.  ]]
    )
    cue = 0
    opt1 =1
    opt2=2
    resp = 1
    w = 0.0
    tau= 1.0
    prob = cfit.prob_choice(cue, opt1, opt2, resp,sr,transition_matrix, w, tau)
    
#     breakpoint()