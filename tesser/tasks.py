#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This function is to evaluate performance in the temporal structure learning, grouping, and inductive inference tasks
- temporal structure learning task
    temp_node_info()

- inductive inference task
    adjacency_mat(node_df)

- grouping task
    path_length(adj_matrix)
"""

import numpy as np
import scipy.spatial.distance as sd
import pandas as pd
from . import util


###########################
#inductive inference task
###########################
#make a function that gets the average inference task data for a participant
def induct_avg(induct_df):
    part_num = int(induct_df["SubjNum"].mean())
    overall_avg = induct_df["Acc"].mean()
    central_trials = induct_df[(induct_df.QuestType == 'Prim')]
    central_avg = central_trials["Acc"].mean()
    bound1_trials = induct_df[(induct_df.QuestType == 'Bound1')]
    bound1_avg = bound1_trials["Acc"].mean()
    bound2_trials = induct_df[(induct_df.QuestType == 'Bound2')]
    bound2_avg = bound2_trials["Acc"].mean()
    return part_num, overall_avg, central_avg, bound1_avg, bound2_avg


#make a function that gets the average inference task data for all participants
def induct_avg_all(data_dir):
    #want to load a list of the participants and loop through
    part_avg_list = []
    part_list = util.subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        this_induct = util.load_induct_df_all(data_dir, part_num)
        this_avg = induct_avg(this_induct)
        part_avg_list.append(this_avg)
    #convert list to df
    part_avg_df = pd.DataFrame(part_avg_list, columns=['Participant', 'Overall', 'Prim', 'Bound1', 'Bound1'])
    return part_avg_df


#make a function that returns a median split of the data based on overall induction performance
def induct_avg_split_high(part_avg_df):
    induct_avg = part_avg_df["Overall"].mean()

    induct_above = part_avg_df[part_avg_df["Overall"] > induct_avg]
    return induct_above


#make a function that returns a median split of the data based on overall induction performance
def induct_avg_split_low(part_avg_df):
    induct_avg = part_avg_df["Overall"].mean()

    induct_below = part_avg_df[part_avg_df["Overall"] < induct_avg]
    return induct_below


###########################
#grouping task
###########################
def group_dist(group_mat, distance='Euclidean'):
    """Calculate object distance from grouping task."""

    rows, cols = np.nonzero(group_mat)
    objects = group_mat[rows, cols]
    coord = np.vstack((rows, cols)).T
    sort_ind = np.argsort(objects)
    coord_sort = coord[sort_ind]
    group_dist_mat = sd.squareform(sd.pdist(coord_sort, distance))
    return group_dist_mat

#make a figure of the grouping task with corresponding colors

#not sure if this function is useful?
def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""

    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat
