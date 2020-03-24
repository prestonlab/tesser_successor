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


import pandas as pd
import numpy as np
import scipy.spatial.distance as sd
import matplotlib
import matplotlib.pyplot as plt
from . import util
from . import network

###########################
# inductive inference task
###########################
# make a function that gets the average inference task data for a participant


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


# make a function that gets the average inference task data for all participants
def induct_avg_all(data_dir):
    # want to load a list of the participants and loop through
    part_avg_list = []
    part_list = util.subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        this_induct = util.load_induct_df_all(data_dir, part_num)
        this_avg = induct_avg(this_induct)
        part_avg_list.append(this_avg)
    # convert list to df
    part_avg_df = pd.DataFrame(part_avg_list, columns=['Participant', 'Overall', 'Prim', 'Bound1', 'Bound1'])
    return part_avg_df


# make a function that returns a median split of the data based on overall induction performance
def induct_avg_split_high(participant_induct_avg_df):
    induct_average = participant_induct_avg_df["Overall"].mean()

    induct_above = participant_induct_avg_df[participant_induct_avg_df["Overall"] > induct_average]
    return induct_above


# make a function that returns a median split of the data based on overall induction performance
def induct_avg_split_low(participant_induct_avg_df):
    induct_average = participant_induct_avg_df["Overall"].mean()

    induct_below = participant_induct_avg_df[participant_induct_avg_df["Overall"] < induct_average]
    return induct_below


###########################
# grouping task
###########################

def group_dist_mat(group_mat):
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
    group_dist_matrix = sd.squareform(sd.pdist(coord_list_sort, 'Euclidean'))
    return group_dist_matrix


# %%

def group_dist_df(group_dist_matrix):
    # first convert to a dataframe
    group_dist_dataframe = pd.DataFrame(group_dist_matrix)

    # reindex the data rows
    group_dist_dataframe.index = group_dist_dataframe.index + 1

    # reindx the data cols
    group_dist_dataframe.columns = group_dist_dataframe.columns + 1

    return group_dist_df  # now index rows, cols = 1:21


# %%

# make unique combo list
def obj_combo(list_objs):
    # have now gotten all the ones in community 1, reindex the data rows:
    list_objs = list_objs.reset_index()
    pair_list = []
    for x in range(len(list_objs)):

        if x == len(list_objs):
            break
        else:
            for q in range(x, len(list_objs) - 1):
                y = q + 1
                pair = [list_objs['node'][x], list_objs['node'][y]]
                pair_list.append(pair)
    return pair_list


# within-community distance between objects in particular community
def specific_within_comm_dist(group_dist_m, pair_list):
    this_dist_tog = []
    for l in range(len(pair_list) - 1):
        this_pair = pair_list[l]
        mat_x = this_pair[0]
        mat_y = this_pair[1]
        # now match to the matrix of similarity
        dist_df = group_dist_df(group_dist_m)
        this_dist = dist_df.loc[mat_x, mat_y]
        this_dist_tog.append(this_dist)

    this_comm_within_dist = np.mean(this_dist_tog)
    return this_comm_within_dist


# %%
# within-community distance between objects
def within_comm_dist(group_dist_matrix):
    within_dist_tog = []
    # getting the structure information
    comm_all_objs = network.temp_node_info()
    # print(structure_info)
    for n in range(1, 4):
        this_comm_objs = comm_all_objs[comm_all_objs["comm"] == n]
        this_comm_pairs = obj_combo(this_comm_objs)
        this_comm_within_dist = specific_within_comm_dist(group_dist_matrix, this_comm_pairs)
        within_dist_tog.append(this_comm_within_dist)
    within_comm_d = np.mean(within_dist_tog)
    return within_comm_d


# %%
def across_comm_dist(group_dist_matrix):
    # getting the structure information
    comm_all_objs = network.temp_node_info()

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
        dist_df = group_dist_df(group_dist_matrix)
        this_dist = dist_df.loc[mat_x, mat_y]
        this_dist_tog.append(this_dist)

    across_comm_distance = np.mean(this_dist_tog)
    return across_comm_distance

def comm_color_map():
    # Color value constants
    colors = {"d_purple": '#7e1e9c',
              "l_purple": '#bf77f6',
              "d_green": '#15b01a',
              "l_green": '#96f97b',
              "d_red": '#e50000',
              "l_red": '#ff474c',
              "grey": '#d8dcd6'}
    
    color_list = []
    for val in range(0, 22):
        # Map the value to a color
        if val == 0:
            color = colors["grey"]
        elif val in comm1_cent_objs['node']:
            color = colors["d_purple"]
        elif val in comm1_bound_objs['node']:
            color = colors["l_purple"]
        elif val in comm2_cent_objs['node']:
            color = colors["d_red"]
        elif val in comm2_bound_objs['node']:
            color = colors["l_red"]
        elif val in comm3_cent_objs['node']:
            color = colors["d_green"]
        elif val in comm3_bound_objs['node']:
            color = colors["l_green"]
        color_list.append(color)       
    return color_list

def plot_dist(group_mat_input):
    color_list = comm_color_map()
    cmap = matplotlib.colors.ListedColormap(color_list)
    plt.imshow(group_mat_input, cmap=cmap)
    
def total_group_num(group_mat_input):
    #make group_mat into 0 and 1 values instead of object numbers
    bin_group_mat = (group_mat > 0).astype(int)
    from skimage import measure
    groupmat_labeled = measure.label(bin_group_mat, connectivity=2, return_num=True)
    labeled_total = groupmat_labeled[1] #total number of groups
    return labeled_total

def group_clusters(group_mat_input):
    #make group_mat into 0 and 1 values instead of object numbers
    bin_group_mat = (group_mat > 0).astype(int)   
    from skimage import measure
    groupmat_labeled = measure.label(bin_group_mat, connectivity=2, return_num=True)
    labeled_array = groupmat_labeled[0] #array of labeled groups
    group_stacked = []
    
    for t in range(1, labeled_total+1):
        this_index = np.where(labeled_array==t)
        this_group = group_mat[this_index] 
        this_group_list = this_group.tolist()
        group_stacked.append(this_group_list)
    return group_stacked

