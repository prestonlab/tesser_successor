import numpy as np
import pandas as pd
import make_data

data, dl = make_data.get_data('/home/rodrigo/Dropbox/tesser_successor/Data/')

def get_objnum(num):

	objnum ={}# object sequence number
	trial = {}# trial number
	rt ={}# responce/reaction time
	for z in range(len(dl)):
		objnum[dl[z]]=data[dl[z]][' objnum']
		# trial[dl[z]]=data[dl[z]][' trial']
		# rt[dl[z]]=data[dl[z]][' rt']
	
	return objnum[dl[num]]
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
		
	
def SRclass_nathum_exp1(envstep, gamma, alpha, p_sample=None):


#     if p_sample==None:
#         p_sample= [.5,.5]
		
	num_states = 21

	SR_agent = SR_Matrix(gamma, alpha, p_sample, num_states)
	s = envstep[0]-1
	
	for entry in envstep[1:]: # go through trajectory till the end
			
		s_new = entry-1
		SR_agent.step (s, s_new)
			
		s = s_new


	return SR_agent.M

	