import pandas
import numpy as np


# receives a tab-separated input file and returns the data frame containing its contents
def read_file (directory):
    with open (directory) as file:
        matrix = pandas.read_csv (file, sep='\t')
    return matrix

# gets the objects column from the experimental run which is inputted
def get_objects(directory):
    matrix = read_file (directory)
    obj_sequence = matrix[[' objnum']]    
    return obj_sequence

def get_induction_data (directory):
    matrix = read_file (directory)
    cue_sequence = np.array(matrix[[' CueNum']])
    opt1_sequence = np.array(matrix[[' Opt1Num']])
    opt2_sequence = np.array(matrix[[' Opt2Num']])
    response_sequence = np.array(matrix[[' Resp']])
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence
    