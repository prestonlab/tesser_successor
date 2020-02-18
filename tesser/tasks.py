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

import os
import pandas as pd
import numpy as np
import scipy.spatial.distance as sd
import matplotlib.pyplot as plt
from util import *
from tasks import *
from network import *

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
    part_list = subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        this_induct = load_induct_df_all(data_dir, part_num)
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

def group_dist_mat(data_dir, subj_num):
    group_mat = load_group(data_dir, subj_num)
    rows, cols = np.nonzero(group_mat)

    # getting the identity of the objects from their coordinates:
    objects = group_mat[rows, cols]

    # stacking and transposing the coordinates:
    # vertically stacking horizontal list of rows, columns, such that row1 = rows, row2 = cols
    coord_stack = np.vstack((rows, cols))
    # then transposing them to two columns, such that col1 = rows, col2 = cols
    coord_list = coord_stack.T
    # recall that the objects list is in the order in which they were interpreted
    # to understand at what index the lowest (1) and highest (21) objects are in, can use np.argsort
    # this gives you the index of each of those objects from low to high as in the objects + coord list
    sort_ind = np.argsort(objects)

    # now can use the coord list and sort by the low to high sort_ind of object index
    # so now the list is the coordinates of objects 1-21 in that order
    coord_list_sort = coord_list[sort_ind]

    # now we can take coord_list_sort (i.e. coordinates of objects in object order 1-21) and make a distance matrix
    group_dist_mat = sd.squareform(sd.pdist(coord_list_sort, 'Euclidean'))
    return group_dist_mat


# %%

def group_dist_df(group_dist_mat):
    # first convert to a dataframe
    group_dist_df = pd.DataFrame(group_dist_mat)

    # reindex the data rows
    group_dist_df.index = group_dist_df.index + 1

    # reindx the data cols
    group_dist_df.columns = group_dist_df.columns + 1

    return group_dist_df  # now index rows, cols = 1:21


# %%

# make unique combo list
def obj_combo(comm_objs):
    # have now gotten all the ones in community 1, reindex the data rows:
    comm_objs = comm_objs.reset_index()
    pair_list = []
    this_dist_tog = []
    for x in range(len(comm_objs)):

        if x == len(comm_objs):
            break
        else:
            for q in range(x, len(comm_objs) - 1):
                y = q + 1
                # print(x, y)
                # print(comm_objs['node'][x], comm_objs['node'][y])
                pair = [comm_objs['node'][x], comm_objs['node'][y]]
                pair_list.append(pair)
    return pair_list


# within-community distance between objects in particular community
def specific_within_comm_dist(comm_objs, data_dir, subj_num):
    # have now gotten all the ones in community 1, reindex the data rows:
    pair_list = obj_combo(comm_objs)
    this_dist_tog = []
    for l in range(len(pair_list) - 1):
        this_pair = pair_list[l]
        mat_x = this_pair[0]
        mat_y = this_pair[1]
        # now match to the matrix of similarity
        dist_mat = group_dist_mat(data_dir, subj_num)
        dist_df = group_dist_df(dist_mat)
        this_dist = dist_df.loc[mat_x, mat_y]
        this_dist_tog.append(this_dist)

    this_comm_within_dist = np.mean(this_dist_tog)
    return (this_comm_within_dist)


# %%
# within-community distance between objects
def within_comm_dist(data_dir, subj_num):
    within_dist_tog = []
    # getting the structure information
    comm_all_objs = temp_node_info()
    # print(structure_info)
    for n in range(1, 4):
        comm_objs = comm_all_objs[comm_all_objs["comm"] == n]
        this_comm_within_dist = specific_within_comm_dist(comm_objs, data_dir, subj_num)
        within_dist_tog.append(this_comm_within_dist)
    within_comm_dist = np.mean(within_dist_tog)
    return (within_comm_dist)


# %%
def across_comm_dist(data_dir, subj_num):
    # getting the structure information
    comm_all_objs = temp_node_info()

    comm1_objs = comm_all_objs[comm_all_objs["comm"] == 1]
    comm2_objs = comm_all_objs[comm_all_objs["comm"] == 2]
    comm3_objs = comm_all_objs[comm_all_objs["comm"] == 3]

    comm1_pairs = obj_combo(comm1_objs)
    comm2_pairs = obj_combo(comm2_objs)
    comm3_pairs = obj_combo(comm3_objs)

    within_comm_pairs = comm1_pairs + comm2_pairs + comm3_pairs
    within_comm_pairs_array = np.array(within_comm_pairs)

    all_comm_pairs = obj_combo(comm_all_objs)
    all_comm_pairs_array = np.array(all_comm_pairs)

    diff = set(map(tuple, all_comm_pairs_array)) - set(map(tuple, within_comm_pairs_array))
    diff_list = list(map(list, diff))

    pair_list = diff_list
    this_dist_tog = []
    for l in range(len(pair_list) - 1):
        this_pair = pair_list[l]
        mat_x = this_pair[0]
        mat_y = this_pair[1]
        # now match to the matrix of similarity
        dist_mat = group_dist_mat(data_dir, subj_num)
        dist_df = group_dist_df(dist_mat)
        this_dist = dist_df.loc[mat_x, mat_y]
        this_dist_tog.append(this_dist)

    across_comm_dist = np.mean(this_dist_tog)
    return across_comm_dist


#make a figure of the grouping task with corresponding colors

#not sure if this function is useful?
def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""

    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat
