import numpy as np
import matplotlib.pyplot as plt
import random
import math

class SR_Matrix():

	''' This class defines a reinforcement learning agent that 
	learns the state-state successor representation without taking actions. 
	Thus, the resulting SR matrix is in the service of prediction. 

	Initalization parameters

	gamma: discount param
	alpha: learning rate
	p_sample: probability of sampling different options, only relevant for testing poilcy dependence
	NUM_STATES: the number of states in the environment to intialize matrices

	Ida Momennejad, 2019'''

	def __init__(self, gamma, alpha, p_sample, NUM_STATES):
		self.gamma = gamma # discount factor
		self.alpha = alpha # learning rate
		self.p_sample = p_sample # p(sampling options)		
		self.M= np.zeros([NUM_STATES, NUM_STATES]) # M: state-state SR    	
		self.onehot=np.eye(NUM_STATES) # onehot matrix, for updating M
		
	def step(self, s, s_new):

		self.update_SR(s, s_new)
        


	def update_SR(self, s, s_new):
		self.M[s] = (1-self.alpha)* self.M[s] + self.alpha * ( self.onehot[s] + self.gamma * self.M[s_new]  )
		
	
        
        
def run_experiment(envstep, gamma, alpha, p_sample=None, verbose=0):

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
        
    num_states = 21

    SR_agent = SR_Matrix(gamma, alpha, p_sample, num_states)
    s = envstep[0]-1
    
    for entry in envstep[1:]: # go through trajectory till the end
            
        s_new = entry-1
        SR_agent.step (s, s_new)
            
        s = s_new


    return SR_agent.M
