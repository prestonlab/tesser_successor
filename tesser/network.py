#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This network function is used to showcase the temporal community structure network used in Tesser
- 21-node temporal community structure in a dataframe
    temp_node_info()

- adjacency matrix for the temporal community structure (i.e. connected nodes)
    adjacency_mat(node_df)

- path legnth between the connected nodes
    path_length(adj_matrix)

- community matrix for the temporal community structure (i.e. assuming just unconnected groups)
    community_mat(node_df)
"""
import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path


def temp_node_info():
    """ INPUT [none]
        OUTPUT node dataframe detailing the community type (1, 2, 3), node type (boundary = 1, non-boundary = 0), and connectedness (1 = connected, 0 = not) of the 21 nodes in temporal community structure  """

    n_node = 21

    # node number
    node = np.arange(1, n_node + 1, dtype=int)
    df = pd.DataFrame(index=node)
    df['node'] = node

    # communities
    df['comm'] = 0
    comm = {1: [1, 2, 3, 18, 19, 20, 21],
            2: [4, 5, 6, 7, 8, 9, 10],
            3: [11, 12, 13, 14, 15, 16, 17]}
    for n, items in comm.items():
        df.loc[items, 'comm'] = n

    # nodetype (0: primary, 1: boundary)
    df['nodetype'] = 0
    df.loc[[3, 18, 4, 10, 11, 17], 'nodetype'] = 1

    # connected community
    df['connect'] = 0
    df.loc[[4, 17], 'connect'] = 1
    df.loc[[3, 11], 'connect'] = 2
    df.loc[[10, 18], 'connect'] = 3
    return df


def community_mat(node_df):
    """ INPUT node_df
        OUTPUT commmunity matrix (i.e. isolated, unconnected) temporal community structure  """

    n_node = node_df.shape[0]
    comm = np.zeros((n_node, n_node), dtype=int)
    for node in node_df.index:
        comm_row = np.zeros(n_node)
        if node_df.loc[node, 'nodetype'] == 0:
            # connected to all items in the community
            comm_row[node_df.comm == node_df.loc[node, 'comm']] = 1
        else:
            # connected to all community nodes except the other boundary
            samecomm = node_df.comm == node_df.loc[node, 'comm']
            isprim = node_df.nodetype == 0
            comm_row[samecomm & isprim] = 1
        comm_row[node - 1] = 0
        comm[node - 1, :] = comm_row
        comm[:, node - 1] = comm_row
    return comm


def adjacency_mat(node_df):
    """ INPUT node_df
        OUTPUT adjacency matrix of temporal community structure  """

    n_node = node_df.shape[0]
    adj = np.zeros((n_node, n_node), dtype=int)
    for node in node_df.index:
        adj_row = np.zeros(n_node)
        if node_df.loc[node, 'nodetype'] == 0:
            # connected to all items in the community
            adj_row[node_df.comm == node_df.loc[node, 'comm']] = 1
        else:
            # connected to boundary node of connected community
            outcon = node_df.comm == node_df.loc[node, 'connect']
            incon = node_df.connect == node_df.loc[node, 'comm']
            adj_row[outcon & incon] = 1

            # connected to all community nodes except the other boundary
            samecomm = node_df.comm == node_df.loc[node, 'comm']
            isprim = node_df.nodetype == 0
            adj_row[samecomm & isprim] = 1
        adj_row[node - 1] = 0
        adj[node - 1, :] = adj_row
        adj[:, node - 1] = adj_row
    return adj


def path_length(adj_mat):
    """ INPUT adjacency matrix of temporal community structure
        OUTPUT shortest path length  """

    n_node = adj_mat.shape[0]
    G = nx.from_numpy_matrix(adj_mat)
    gpath = shortest_path(G)
    pathlen = np.zeros((n_node, n_node))
    for source, pathdict in gpath.items():
        for target, pathlist in pathdict.items():
            pathlen[source, target] = len(pathlist) - 1
    return pathlen
