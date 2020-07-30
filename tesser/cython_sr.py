#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Learning module for tesser simulations. Functions for learning/running
experiment phases are based on Ida Momennejad's state-state successor
representation learning agent.

These include:

- Class that defines a reinforcement learning agent which learns the state-state
  successor representation without taking actions.
    SRMatrix()
- Function uses the reinforcement learning agent class in SRMatrix to learn.
    run_experiment(envstep, gamma, alpha, M, n_states)

- Train an SR matrix on the structure-learning task.
    learn_sr(df, gamma, alpha)

- Computes the matrix to which SR learning should converge, by summing a
    geometric matrix series.
    compute_limit_matrix(gamma, adjacency, n_states)

- Computes the correlation matrix for a matrix's rows.
    correlate_rows(matrix)

- Computes the correlation matrix for a matrix's columns.
    correlate_columns(matrix)

- Computes the norm or correlation between the SR matrix and the limit matrix,
  given a subject & values for gamma, alpha
    compute_correlations(df, option, gamma, alpha)

"""
import numpy as np
from . import csr

def learn_sr(df, gamma, alpha, n_states):
    """
    Train an SR matrix on the structure-learning task.

    Parameters
    ----------
    df : pandas.DataFrame
        Structure learning task trials. Must have fields:
        objnum - object number (starting from 1)
        part - part number
        run - run number

    gamma : float
        Discounting factor.

    alpha : float
        Learning rate.
        
    n_states : int
        The number of states in the environment to initialize matrices.

    Returns
    -------
    M : numpy.array
        SR Matrix for all parts and run for a given subject.
    """
    M = np.zeros([n_states, n_states])
    onehot= np.eye(n_states, dtype = np.dtype('i'))

    envstep = df.objnum.to_numpy() -1
    envstep = envstep.astype(np.dtype('i'))
    csr.SR(envstep, gamma, alpha, M, n_states, onehot)
    return M
