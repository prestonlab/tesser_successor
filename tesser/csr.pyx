cimport cython


@cython.boundscheck(False)
@cython.wraparound(False) 
def SR(int [:] envstep, double gamma, double alpha, double [:,:] M, int n_states, int [:,:] onehot):
    cdef int s, s_new, i, n
    cdef Py_ssize_t n_steps = envstep.shape[0]
    # set initial state
    s = envstep[0]
    for n in range(1,n_steps):
        s_new = envstep[n]
        # update matrix based on state transition
        # SR_agent.step(s, s_new)
        for i in range(n_states):
            M[s,i] = (1 - alpha) * M[s, i] + alpha * (
                onehot[s_new, i] + gamma * M[s_new, i]
            )
        s = s_new
