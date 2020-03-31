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
from scipy.spatial.distance import pdist, squareform
import matplotlib
import matplotlib.pyplot as plt
from . import util
from . import network


###########################
# inductive inference task
###########################


def induct_avg(induct_df):
    """Get average induction task performance by question type."""

    part_num = int(induct_df["SubjNum"].mean())
    overall_avg = induct_df["Acc"].mean()
    central_trials = induct_df[(induct_df.QuestType == 'Prim')]
    central_avg = central_trials["Acc"].mean()
    bound1_trials = induct_df[(induct_df.QuestType == 'Bound1')]
    bound1_avg = bound1_trials["Acc"].mean()
    bound2_trials = induct_df[(induct_df.QuestType == 'Bound2')]
    bound2_avg = bound2_trials["Acc"].mean()
    return part_num, overall_avg, central_avg, bound1_avg, bound2_avg


def induct_avg_all(data_dir):
    """Average induction task performance for all participants."""

    # want to load a list of the participants and loop through
    part_avg_list = []
    part_list = util.subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        this_induct = util.load_induct_subject(data_dir, part_num)
        this_avg = induct_avg(this_induct)
        part_avg_list.append(this_avg)
    # convert list to df
    part_avg_df = pd.DataFrame(part_avg_list, columns=['participant', 'overall', 'prim', 'bound1', 'bound2'])
    return part_avg_df


def induct_avg_split_high(all_participant_induct_avg_df):
    """High performers based on a median split of induction."""

    induct_average = all_participant_induct_avg_df["overall"].mean()

    induct_above = all_participant_induct_avg_df[all_participant_induct_avg_df["overall"] > induct_average]
    return induct_above


def induct_avg_split_low(all_participant_induct_avg_df):
    """Low performers based on a median split of induction."""

    induct_average = all_participant_induct_avg_df["overall"].mean()

    induct_below = all_participant_induct_avg_df[all_participant_induct_avg_df["overall"] < induct_average]
    return induct_below


def induct_bias(induct_df):
    """Get average induction task bias by question type."""
    full_len = len(induct_df)
    #primary
    induct_prim = induct_df[induct_df["QuestType"]=="Prim"]
    prim_len = len(induct_prim)
    #bound1
    induct_bound1 = induct_df[induct_df["QuestType"]=="Bound1"]
    bound1_len = len(induct_bound1)
    #bound2
    induct_bound2 = induct_df[induct_df["QuestType"]=="Bound2"]
    bound2_len = len(induct_bound2)

    #get overall mean
    overall_acc_total = len(induct_df[induct_df["Acc"]==1])
    overall_acc_prop = overall_acc_total/full_len
    overall_mis_total = len(induct_df[induct_df["Acc"]==0])
    overall_mis_prop = overall_mis_total/full_len

    overall_bias = overall_acc_prop - overall_mis_prop
    
    #get primary mean (acc + RT)
    prim_acc_total = len(induct_prim[induct_prim["Acc"]==1])
    prim_acc_prop = prim_acc_total/prim_len
    prim_mis_total = len(induct_prim[induct_prim["Acc"]==0])
    prim_mis_prop = prim_mis_total/prim_len

    prim_bias = prim_acc_prop - prim_mis_prop

    #get bound1 mean (acc + RT)
    bound1_acc_total = len(induct_bound1[induct_bound1["Acc"]==1])
    bound1_acc_prop = bound1_acc_total/bound1_len
    bound1_mis_total = len(induct_bound1[induct_bound1["Acc"]==0])
    bound1_mis_prop = bound1_mis_total/bound1_len

    bound1_bias = bound1_acc_prop - bound1_mis_prop  
    
    #get bound2 mean (acc + RT)
    bound2_acc_total = len(induct_bound2[induct_bound2["Acc"]==1])
    bound2_acc_prop = bound2_acc_total/bound2_len
    bound2_mis_total = len(induct_bound2[induct_bound2["Acc"]==0])
    bound2_mis_prop = bound2_mis_total/bound2_len

    bound2_bias = bound2_acc_prop - bound2_mis_prop
    
    return overall_bias, prim_bias, bound1_bias, bound2_bias


def induct_bias_all(data_dir):
    """Average induction task bias for all participants."""

    # want to load a list of the participants and loop through
    part_bias_list = []
    part_list = util.subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        this_induct = util.load_induct_subject(data_dir, part_num)
        overall_bias, prim_bias, bound1_bias, bound2_bias = induct_bias(this_induct)
        part_bias = [part_num, overall_bias, prim_bias, bound1_bias, bound2_bias]        
        part_bias_list.append(part_bias)
    # convert list to df
    part_bias_df = pd.DataFrame(part_bias_list, columns=['participant', 'overall', 'prim', 'bound1', 'bound2'])
    return part_bias_df

###########################
# grouping task
###########################

def euc_dist_mat(data_dir, subject):
    """Calculate pairwise object distances from a grouping array output."""
    group_load = util.load_group_subject(data_dir, subject)
    df_group = pd.DataFrame(group_load)
    df_group.columns = np.arange(len(df_group.columns))
    
    all_num = []
    all_row = []
    all_col = []
    for n in range(1, 22, 1):
        this_row_col = np.where(df_group == n) #row [0], column [1]
        this_row = np.int(this_row_col[0])
        this_col = np.int(this_row_col[1])
        all_num.append(n)
        all_row.append(this_row)
        all_col.append(this_col)
  
    all_num_row_col = pd.DataFrame(list(zip(all_num, all_row, all_col)), columns =['node', 'row', 'col'])
    all_row_col = pd.DataFrame(list(zip(all_row, all_col)), columns =['row', 'col'])
    all_row_col_array = all_row_col.values
    
    dist_bw_sq = squareform(pdist(all_row_col_array, "euclidean"))
    return dist_bw_sq


def euc_dist_comm_avg(data_dir, subject):
    # load structure learning and induction data for one subject
    dist_bw_sq = euc_dist_mat(data_dir, subject)

    #get only the upper half, excluding diagonal (b/c the same above and below diagonal)
    dist_upper = np.triu(dist_bw_sq)

    #add in NaNs next
    dist_upper[np.tril_indices(dist_upper.shape[0], 0)] = np.nan

    comm_struct = network.temp_node_info()

    #now create an adjacency matrix basd on the relationships b/w items
    #can use this adjacency matrix to create a mask to cover the distance above 
    comm = network.community_mat(comm_struct)
    
    #get only the upper half, excluding diagonal (b/c the same above and below diagonal)
    comm_upper = np.triu(comm)

    #add in 7 for the diagonal that does not matter, b/c rest will be using 1 and 0 for (within and across) 
    comm_upper[np.tril_indices(comm_upper.shape[0], 0)] = .01

    #now to calculate the distances whose adj_upper == 1 (within) and adj_upper == 0 (across)
    #within
    comm_within_dist = np.ma.masked_where(comm_upper == 0, dist_upper)

    #now to take the mean of the distances
    within_avg_dist_comm = np.nanmean(comm_within_dist)

    #across
    comm_across_dist = np.ma.masked_where(comm_upper == 1, dist_upper)

    #now to take the mean of the distances
    across_avg_dist_comm = np.nanmean(comm_across_dist)

    return within_avg_dist_comm, across_avg_dist_comm


def euc_dist_comm_avg_all(data_dir):
    """Average induction task performance for all participants."""

    # want to load a list of the participants and loop through
    part_euc_avg_list = []
    part_list = util.subj_list()
    for i in range(len(part_list)):
        part_num = part_list[i]
        within_dist, across_dist = euc_dist_comm_avg(data_dir, part_num)
        part_avg = [part_num, within_dist, across_dist]
        part_euc_avg_list.append(part_avg)
    # convert list to df
    part_euc_avg_df = pd.DataFrame(part_euc_avg_list, columns=['participant', 'within_dist', 'across_dist'])
    return part_euc_avg_df


def comm_color_map():
    
    #community members
    comm_all_objs = network.temp_node_info()
    comm1_cent_objs = comm_all_objs[(comm_all_objs["comm"] == 1) & (comm_all_objs["nodetype"] == 0)]
    comm1_bound_objs = comm_all_objs[(comm_all_objs["comm"] == 1) & (comm_all_objs["nodetype"] == 1)]
    comm2_cent_objs = comm_all_objs[(comm_all_objs["comm"] == 2) & (comm_all_objs["nodetype"] == 0)]
    comm2_bound_objs = comm_all_objs[(comm_all_objs["comm"] == 2) & (comm_all_objs["nodetype"] == 1)]
    comm3_cent_objs = comm_all_objs[(comm_all_objs["comm"] == 3) & (comm_all_objs["nodetype"] == 0)]
    comm3_bound_objs = comm_all_objs[(comm_all_objs["comm"] == 3) & (comm_all_objs["nodetype"] == 1)]

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


def plot_dist(group_array):
    color_list = comm_color_map()
    cmap = matplotlib.colors.ListedColormap(color_list)
    plt.imshow(group_array, cmap=cmap)


def total_group_num(group_array):
    # make group_mat into 0 and 1 values instead of object numbers
    bin_group_mat = (group_array > 0).astype(int)
    from skimage import measure
    groupmat_labeled = measure.label(bin_group_mat, connectivity=2, return_num=True)
    labeled_total = groupmat_labeled[1] #total number of groups
    return labeled_total


def group_clusters(group_array):
    # make group_mat into 0 and 1 values instead of object numbers
    bin_group_mat = (group_array > 0).astype(int)   
    from skimage import measure
    groupmat_labeled = measure.label(bin_group_mat, connectivity=2, return_num=True)
    labeled_total = total_group_num(group_array)
    labeled_array = groupmat_labeled[0] #array of labeled groups
    group_stacked = []
    
    for t in range(1, labeled_total+1):
        this_index = np.where(labeled_array==t)
        this_group = group_array[this_index] 
        this_group_list = this_group.tolist()
        group_stacked.append(this_group_list)
    return group_stacked
