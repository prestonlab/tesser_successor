import numpy as np
import pandas as pd
import make_data

data, dl, induc_data, idl = make_data.get_data('/home/rodrigo/Dropbox/tesser_successor/Data/')

def get_objnum(num):

    objnum ={}# object sequence number
    trial = {}# trial number
    rt ={}# responce/reaction time
    for z in range(len(dl)):
        objnum[dl[z]]=data[dl[z]][' objnum']
        # trial[dl[z]]=data[dl[z]][' trial']
        # rt[dl[z]]=data[dl[z]][' rt']

    return objnum[dl[num]]

def get_induc_data(num):

    CueNum = {}# cue sequence number
    Opt1Num = {}# object sequence number for option 1
    Opt2Num = {}# object sequence number for option 2
    Resp = {} # response number

    for z in range(len(idl)):

        CueNum[idl[z]]=induc_data[idl[z]][' CueNum']
        Opt1Num[idl[z]]=induc_data[idl[z]][' Opt1Num']
        Opt2Num[idl[z]]=induc_data[idl[z]][' Opt2Num']
        Resp[idl[z]]=induc_data[idl[z]][' Resp']
        
    return CueNum[idl[num]], Opt1Num[idl[num]], Opt2Num[idl[num]], Resp[idl[num]]
        
    # Derived from explore_simulations.py by Demetrius Rowland
def make_envstep(direct_num):
    envstep = []
    objects = get_objnum(direct_num)
    index = 0

    for index in range ( len(objects)):
        obj = int(objects[index])
        try:
            envstep.append (obj)
        except:
            pass
    return envstep


'''
 SR_no_action reduced version
'''
reward = 1 # for our purposes we will keep the reward parameter constant

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

    def __init__(self, gamma, alpha, NUM_STATES, M):
        self.gamma = gamma # discount factor
        self.alpha = alpha # learning rate
        self.M= np.zeros([NUM_STATES, NUM_STATES]) # M: state-state SR    	
        self.onehot=np.eye(NUM_STATES) # onehot matrix, for updating M

    def step(self, s, s_new):

        self.update_SR(s, s_new)



    def update_SR(self, s, s_new):
        self.M[s] = (1-self.alpha)* self.M[s] + self.alpha * ( self.onehot[s] + self.gamma * self.M[s_new]  )


def SRclass_nathum_exp1(envstep, gamma, alpha, M):

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

    num_states = 21

    SR_agent = SR_Matrix(gamma, alpha, num_states, M)
    s = envstep[0]-1

    for entry in envstep[1:]: # go through trajectory till the end

        s_new = entry-1
        SR_agent.step (s, s_new)

        s = s_new


    return SR_agent.M

