import numpy as np
import matplotlib.pyplot as plt
import random
from SR_no_action import SR_no_action
#from dyna_replay import dyna_replay

def SRclass_nathum_exp1(envstep, gamma, alpha, p_sample=None, verbose=0):

    ''' This function uses the reinfrocement learning agent class in 
        SR_no_action.py to learn.
        Here the function takes the environment from Experiment 1 in our
        Nat Hum Beh paper & learns predictive representations with the 
        specified learning rate and scale. 

        Note: This is not SR dyna, nor SR-MB. 
        This agent only learns the SR.

        Inputs: 
        
        envstep:  generated with ida_envs.generate_nathumbeh_env1()        
        gamma: discount parameter, determines scale of predictive representations
        alpha: learning rate
        p_sample: prob sampling each of the two sequences 
        verbose: TRUE or FALSE 

        Outputs:

        M: SR matrix 
        W: value weights W
        memory: memory of episodes 
        episodies: # episodes it takes to reach convergence 

        Ida Momennejad, NYC, 2019'''

    if p_sample==None:
        p_sample= [.5,.5]
        
    num_states = 22

    SR_agent = SR_no_action(gamma, alpha, p_sample, num_states)
    done = False
    
    for time in range (1000):
        # sample a starting point [really determined by experiment]
        s = envstep[0]

        done = False
        index = 0
        while not done: # go through trajectory till the end
            
            s_new, done = envstep [index]
            SR_agent.step (s, s_new)
            
            s = s_new
            index += 1


    return SR_agent.M


