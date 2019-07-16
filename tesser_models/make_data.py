import numpy as np
import pandas as pd
import os


def get_data(path):
    
    # simple for loop to obtain all data .txt files including subdirectories.
    files = []
    # r=root, d=directories, f = files
    for r, d, f in os.walk(path):
        for file in f:
            if '.txt' in file:
                files.append(os.path.join(r, file))

    # list of file paths for the required data            
    listd=[]
    for f in files:
        listd.append(f)
    listd.sort()

    # list of data path labels to use as keys
    dl=[]
    # dictionary of data sets with data labels
    data={}


    for i in range (len(listd)):
        s = listd[i]
        lb = s[44:]
        dl.append(lb)
        Test = pd.read_csv(s,sep="\t")#[' objnum']
        Test.replace([" NaN"], np.nan, inplace = True) # Drops NaN values from DataFrame
        Test = Test.dropna()
        b = Test.reset_index(drop=True) # Resets the index to start at 0
        data[lb] = b

    return data, dl
