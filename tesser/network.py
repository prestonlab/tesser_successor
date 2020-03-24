#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Information about the temporal community structure used in Tesser.

21-node temporal community structure in a dataframe:
    temp_node_info()

Adjacency matrix for the temporal community structure (i.e. connected
nodes):
    adjacency_mat(node_df)

Path legnth between the connected nodes:
    path_length(adj_matrix)

Community matrix for the temporal community structure (i.e. assuming
just unconnected groups):
    community_mat(node_df)
"""

import numpy as np
import pandas as pd
import networkx as nx
from . import sr
from networkx.algorithms.shortest_paths.generic import shortest_path
import matplotlib.pyplot as plt


def temp_node_info():
    """Get information about nodes in the structured learning task network.

    Returns
    -------
    pandas.DataFrame
        Dataframe detailing the community type (1, 2, 3), node type
        (boundary = 1, non-boundary = 0), and connectedness
        (1 = connected, 0 = not) of the 21 nodes in temporal community
        structure
    """

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
    """Label adjacent node pairs that are in the same community.

    Parameters
    ----------
    node_df : pandas.DataFrame
        Information about node structure. See temp_node_info.

    Returns
    -------
    comm : array
        A nodes x nodes array with 1 for node pairs that are
        adjacent and in the same community and 0 for node pairs that
        are non-adjacent or in different communities.
    """

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
    """Create an adjacency matrix based on network node information.

    Parameters
    ----------
    node_df : pandas.DataFrame
        Information about node structure. See temp_node_info.

    Returns
    -------
    adj : array
        A nodes x nodes array with 1 for adjacent node pairs and 0 for
        unconnected nodes.
    """

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
    """Calculate shortest path length from an adjacency matrix."""

    n_node = adj_mat.shape[0]
    G = nx.from_numpy_matrix(adj_mat)
    gpath = shortest_path(G)
    pathlen = np.zeros((n_node, n_node))
    for source, pathdict in gpath.items():
        for target, pathlist in pathdict.items():
            pathlen[source, target] = len(pathlist) - 1
    return pathlen


def plot_adjacency_matrix():
    nodes = temp_node_info()
    adjacency = adjacency_mat(nodes)
    L = sr.compute_limit_matrix(0.5, adjacency, 21)
    plt.matshow(L, vmin=0, vmax=0.5)
    plt.colorbar()