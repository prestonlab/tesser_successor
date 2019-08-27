# functions to read behavioral data (structure learning, induction)
import numpy as np
import pandas as pd
import os

p = '/home/rodrigo/Dropbox/tesser_successor/Data/'

def read_file(path=p):
    # simple for loop to obtain all data .txt files including subdirectories.
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.txt' in file:
                files.append(os.path.join(r, file))
    # list of file paths for the required data            
    listd = []
    for f in files:
        listd.append(f)
    listd.sort()
    # list of data path labels to use as keys
    induc_list = []
    # dictionary of induction data sets with directory labels
    induc_data = {}  
    # list of data path labels to use as keys
    struc_list = []
    # dictionary of data sets with data labels
    struc_data = {}
    # for loop design to differentiate data between structured learning and generalized induction
    for i in range (len(listd)):
        if listd[i][-len('InductGen.txt'):] == 'InductGen.txt':
            s = listd[i]
            b = s[44:]
            induc_list.append(b)
            matrix = pd.read_csv(s,sep="\t")#[' objnum']
            induc_data[b] = matrix
        else:     
            s = listd[i]
            lb = s[44:]
            struc_list.append(lb)
            matrix = pd.read_csv(s,sep="\t")#[' objnum']
            matrix.replace([" NaN"], np.nan, inplace = True) # Drops NaN values from DataFrame
            matrix = matrix.dropna()
            matrix = matrix.reset_index(drop=True) # Resets the index to start at 0
            struc_data[lb] = matrix
    return struc_data, struc_list, induc_data, induc_list

# input index of desired data with path to the data files returns object sequence for structured learning runs
def get_objnum(num, path = p):
    objnum ={}# object sequence number
    data, dl, dummy2, dummy3 = read_file(path)
    for z in range(len(dl)):
        objnum[dl[z]]=data[dl[z]][' objnum']
    return objnum[dl[num]]

# input directory for desired data with path to the data files returns object sequence for structured learning runs
def get_objects_dir(directory, path = p):
    data, dummy1, dummy2, dummy3 = read_file(path)
    matrix = data[directory]
    obj_sequence = matrix[[' objnum']]    
    return obj_sequence

# input index of desired data with path to the data files returns four variables from the generalized induction data
def get_induction_data(num, path = p):
    CueNum = {}# cue sequence number
    Opt1Num = {}# object sequence number for option 1
    Opt2Num = {}# object sequence number for option 2
    Resp = {} # response number
    dummy1, dummy2, induc_data, idl = read_file(path)
    for z in range(len(idl)):
        CueNum[idl[z]]=induc_data[idl[z]][' CueNum']
        Opt1Num[idl[z]]=induc_data[idl[z]][' Opt1Num']
        Opt2Num[idl[z]]=induc_data[idl[z]][' Opt2Num']
        Resp[idl[z]]=induc_data[idl[z]][' Resp']  
    return CueNum[idl[num]], Opt1Num[idl[num]], Opt2Num[idl[num]], Resp[idl[num]]

# input directory for desired data with path to the data files returns four variables from the generalized induction data
def get_induction_data_dir(directory, path = p):
    dummy1, dummy2, data, dummy3 = read_file(path)
    matrix  = data[directory]
    cue_sequence = np.array(matrix[[' CueNum']])
    opt1_sequence = np.array(matrix[[' Opt1Num']])
    opt2_sequence = np.array(matrix[[' Opt2Num']])
    response_sequence = np.array(matrix[[' Resp']])
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence

# drops rows containing NaN values resets index
def drop_nan(data):
    data.replace([" NaN"], np.nan, inplace = True) # Drops NaN values from DataFrame
    data = data.dropna()
    data = data.reset_index(drop=True) # Resets the index to start at 0
    return data

# Based on subject number, function returns the index of directory where the subject data starts 
def initial_subj_run(num, path = p):
    count = []
    dummy, dl, dummy2, dummy3 = read_file(p)
    for i in range(0, len(dl), 11):
        count.append(i)
    return count[num]

# Print the number of participants in the data set and a list of subject names
def subject_list(path = p):
    dummy1, dummy2, dummy3, idl = read_file(path)
    print("There are %d subjects in this data set"% (len(idl)))
    print()
    for x in range(len(idl)):
        print('Subject %s is named %s'% (x,idl[x][:17]))
    