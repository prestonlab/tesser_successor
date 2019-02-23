"""Module for tesser-related utilities."""

import numpy as np
import pandas as pd
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path

def node_info():
    """Return a dataframe with basic node information."""

    n_node = 21

    # node number
    node = np.arange(1, n_node+1, dtype=int)
    df = pd.DataFrame(index=node)

    # communities
    df['comm'] = 0
    comm = {1:[1,2,3,18,19,20,21],
            2:[4,5,6,7,8,9,10],
            3:[11,12,13,14,15,16,17]}
    for n, items in comm.items():
        df.loc[items,'comm'] = n

    # nodetype (0: primary, 1: boundary)
    df['nodetype'] = 0
    df.loc[[3,18,4,10,11,17],'nodetype'] = 1

    # connected community
    df['connect'] = 0
    df.loc[[4,17],'connect'] = 1
    df.loc[[3,11],'connect'] = 2
    df.loc[[10,18],'connect'] = 3
    return df

def adjacency(df):
    """Determine adjacency matrix from node information."""

    n_node = df.shape[0]
    adj = np.zeros((n_node, n_node), dtype=int)
    for node in df.index:
        adj_row = np.zeros(n_node)
        if df.loc[node,'nodetype'] == 0:
            # connected to all items in the community
            adj_row[df.comm == df.loc[node,'comm']] = 1
        else:
            # connected to boundary node of connected community
            outcon = df.comm == df.loc[node,'connect']
            incon = df.connect == df.loc[node,'comm']
            adj_row[outcon & incon] = 1

            # connected to all community nodes except the other boundary
            samecomm = df.comm == df.loc[node,'comm']
            isprim = df.nodetype == 0
            adj_row[samecomm & isprim] = 1
        adj_row[node-1] = 0
        adj[node-1,:] = adj_row
        adj[:,node-1] = adj_row
    return adj

def path_length(adj):
    """Determine the shortest path between each pair of nodes."""
    
    n_node = adj.shape[0]
    G = nx.from_numpy_matrix(adj)
    gpath = shortest_path(G)
    pathlen = np.zeros((n_node, n_node))
    for source, pathdict in gpath.items():
        for target, pathlist in pathdict.items():
            pathlen[source,target] = len(pathlist) - 1            
    return pathlen
