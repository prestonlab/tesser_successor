#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Learning module for tesser simulations. Functions for learning/running
experiment phases are based on Ida Momennejad's state-state successor
representation learning agent.

These include:

- Train an SR matrix on the structure-learning task.
    learn_sr(df, gamma, alpha)
    
"""
import numpy as np
from . import csr
import os

def learn_sr(df, gamma, alpha, n_states):
    """
    Train an SR matrix on the structure-learning task.

    Parameters
    ----------
    df : pandas.DataFrame
        Structure learning task trials. Must have fields:
        objnum - object number (starting from 1)

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

    if 'objnum' in df:
        envstep = df.objnum.to_numpy() -1
    else:
        envstep = df['object'].to_numpy() - 1
    envstep = envstep.astype(np.dtype('i'))
    csr.SR(envstep, gamma, alpha, M, n_states, onehot)
    return M

def transition_indiv(subj_df, n_states):
    df=subj_df.reset_index()
    df =df.objnum
    matrix = np.zeros([n_states, n_states])
    for i, j in enumerate(df):
        try:
            matrix[j-1,df[i+1]-1]+=1
        except:
            KeyError
    for row in range(n_states):
        matrix[row] /= np.sum(matrix[row])
    return matrix

def transition_all(struct_df, n_states):
    dict_array = {}
    subjects = struct_df.SubjNum.unique()
    for subject in subjects:
        subj_filter = f'SubjNum == {subject}'
        subj_struct = struct_df.query(subj_filter)
        dict_array[subject] = transition_indiv(subj_struct, n_states)
    return dict_array

def transition_multiple(subjects, directory):
    '''
    subjects = list of integers representing the number assigned to each participant
    directory = string that describes the path for where to find trasitional matrices
    
    '''
    
    return  {subject: np.loadtxt(os.path.join(directory, f"Subject_{subject}_tmatrix.txt")) for subject in subjects}