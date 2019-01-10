# get_adj_tesser returns the adjacency matrix for the graph of experiment
def get_a():
    import pandas as pd
    import numpy as np
    from numpy.linalg import inv
    import matplotlib.pyplot as plt
    
    #create an identity matrix, I
    I = np.identity(21)
    #create the adjacency matrix, A, where 1s will populate every connected node in 21-node temporal community structure
    A = np.zeros((21, 21))

    # TesserScan community memberships by node, 
    # last two numbers = boundary nodes
    # subtracting each community membership by 1 
    # since Python array index starts at 0
    comm1_prim = np.array([1, 2, 19, 20, 21]); 
    comm1_prim_py = comm1_prim - 1; 
    comm1_bound = np.array([3, 18])
    comm1_bound_py = comm1_bound - 1;
    connected_bound1 = np.array([3, 4])
    connected_bound1_py = connected_bound1 - 1;

    comm2_prim = np.array([5, 6, 7, 8, 9])
    comm2_prim_py = comm2_prim - 1;
    comm2_bound = np.array([4, 10])
    comm2_bound_py = comm2_bound - 1; 
    connected_bound2 = np.array([10, 11])
    connected_bound2_py = connected_bound2 - 1; 

    comm3_prim = np.array([12, 13, 14, 15, 16])
    comm3_prim_py = comm3_prim - 1
    comm3_bound = np.array([11, 17])
    comm3_bound_py = comm3_bound - 1
    connected_bound3 = np.array([17, 18])
    connected_bound3_py = connected_bound3 - 1

    # connecting primary and boundary nodes of community
    #community1 
    for r in range(comm1_prim_py.size):
        ro = comm1_prim_py[r]
        for c in range(comm1_prim_py.size):
            co = comm1_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm1_bound_py.size):
        ro = comm1_bound_py[r]
        for c in range(comm1_prim_py.size):
            co = comm1_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm1_prim_py.size):
        ro = comm1_prim_py[r]
        for c in range(comm1_bound_py.size):
            co = comm1_bound_py[c]
            A[ro, co] = 1
        
    for r in range(connected_bound1_py.size):
        ro = connected_bound1_py[r]
        for c in range(connected_bound1_py.size):
            co = connected_bound1_py[c]
            A[ro, co] = 1  
        
        
#community2 
    for r in range(comm2_prim_py.size):
        ro = comm2_prim_py[r]
        for c in range(comm2_prim_py.size):
            co = comm2_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm2_bound_py.size):
        ro = comm2_bound_py[r]
        for c in range(comm2_prim_py.size):
            co = comm2_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm2_prim_py.size):
        ro = comm2_prim_py[r]
        for c in range(comm2_bound_py.size):
            co = comm2_bound_py[c]
            A[ro, co] = 1
        
    for r in range(connected_bound2_py.size):
        ro = connected_bound2_py[r]
        for c in range(connected_bound2_py.size):
            co = connected_bound2_py[c]
            A[ro, co] = 1        
        
        
#community3 
    for r in range(comm3_prim_py.size):
        ro = comm3_prim_py[r]
        for c in range(comm3_prim_py.size):
            co = comm3_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm3_bound_py.size):
        ro = comm3_bound_py[r]
        for c in range(comm3_prim_py.size):
            co = comm3_prim_py[c]
            A[ro, co] = 1
        
    for r in range(comm3_prim_py.size):
        ro = comm3_prim_py[r]
        for c in range(comm3_bound_py.size):
            co = comm3_bound_py[c]
            A[ro, co] = 1
        
    for r in range(connected_bound3_py.size):
        ro = connected_bound3_py[r]
        for c in range(connected_bound3_py.size):
            co = connected_bound3_py[c]
            A[ro, co] = 1      
        
    #make sure that all the diagonals (i.e. items to iteslf) are 0
    np.fill_diagonal(A, 0)

    return A
