#!/usr/bin/env python

from tesser import util
from tesser import network 
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy import stats
import matplotlib.pyplot as plt

def euc_dist_mat(data_dir, subject):
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

def euc_dist_adj_avg(data_dir, subject):
    # load structure learning and induction data for one subject
    dist_bw_sq = euc_dist_mat(data_dir, subject)

    #get only the upper half, excluding diagonal (b/c the same above and below diagonal)
    dist_upper = np.triu(dist_bw_sq)

    #add in NaNs next
    dist_upper[np.tril_indices(dist_upper.shape[0], 0)] = np.nan

    comm_struct = network.temp_node_info()

    #now create an adjacency matrix basd on the relationships b/w items
    #can use this adjacency matrix to create a mask to cover the distance above 
    adj = adjacency(comm_struct)

    #get only the upper half, excluding diagonal (b/c the same above and below diagonal)
    adj_upper = np.triu(adj)

    #add in 7 for the diagonal that does not matter, b/c rest will be using 1 and 0 for (within and across) 
    adj_upper[np.tril_indices(adj_upper.shape[0], 0)] = .01

    #now to calculate the distances whose adj_upper == 1 (within) and adj_upper == 0 (across)
    #within
    adj_within_dist = np.ma.masked_where(adj_upper == 0, dist_upper)

    #now to take the mean of the distances
    within_avg_dist_adj = np.nanmean(adj_within_dist)

    #across
    adj_across_dist = np.ma.masked_where(adj_upper == 1, dist_upper)

    #now to take the mean of the distances
    across_avg_dist_adj = np.nanmean(adj_across_dist)

    return within_avg_dist_adj, across_avg_dist_adj

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
