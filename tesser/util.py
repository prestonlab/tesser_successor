#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# functions to read behavioral data (structured learning, induction, grouping data)
import numpy as np
import pandas as pd
import os
from glob import glob

#load subject directory  
def get_subj_dir(data_dir, subject):
    """Get the data directory for a given subject number."""

    # check that the base directory exists
    if not os.path.exists(data_dir):
        raise IOError(f'Directory does not exist: {data_dir}')

    # look for directories with the correct pattern
    dir_search = glob(os.path.join(data_dir, f'tesserScan_{subject}_*'))
    if len(dir_search) != 1:
        raise IOError(f'Problem finding subject directory for {subject}')

    return dir_search[0]

#load structured learning by a subject, part (1: non-scanned, 2: scanned), particular run  
def load_struct_run(data_dir, subject, part, run):
    """Load data for one structured learning run."""

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject)

    # search for a file with the correct name formatting
    file_pattern = f'tesserScan_{subject}_*_StructLearn_Part{part}_Run_{run}.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject}, part {part}, run {run}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = pd.read_csv(run_file, sep='\t', skipinitialspace=True)

    # add a field indicating the experiment part
    df['part'] = part

    # remove trials with no object (these are just filler trials)
    df = df.loc[np.logical_not(np.isnan(df.objnum)), :]

    # convert object labels to integer
    df = df.astype({'objnum': 'int'})

    return df

#load structured learning in entirety  
def load_struct(data_dir, subject):
    """Load all structured learning data for a subject."""

    # list of all runs to load
    parts = (1, 2)
    part_runs = {1: range(1, 6), 2: range(1, 7)}

    # load individual runs
    df_list = []
    for part in parts:
        for run in part_runs[part]:
            run_df = load_struct_run(data_dir, subject, part, run)
            df_list.append(run_df)

    # concatenate into one data frame
    df = pd.concat(df_list, sort=False)
    return df

#load inductive generalization data by subject 
def load_induction(data_dir, subject):
    """Load data generalized induction data by subject."""

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject)

    # search for a file with the correct name formatting
    file_pattern = f'tesserScan_{subject}_*_InductGen.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = pd.read_csv(run_file, sep='\t', skipinitialspace=True)

    return df

#load inductive generalization by trial by trial questions
def get_induction_data(DFRAME):
    """ INPUT DataFrame 
       OUTPUT four variable sequences for generalized induction """
    data = DFRAME
    cue_sequence = data["CueNum"].values
    opt1_sequence = data["Opt1Num"].values
    opt2_sequence = data["Opt2Num"].values
    response_sequence = data["Resp"].values
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence

#drop the NaN values from inductive generalization questions 
def drop_nan(DATA):
    """  Drops NaN values from DataFrame """
    DATA.replace(["NaN"], np.nan, inplace=True)
    DATA = DATA.dropna()
    DATA = DATA.reset_index(drop=True)  # Resets the index to start at 0
    return DATA

#load grouping data by subject 
def load_group(data_dir, subject):
    """Load data generalized induction data by subject."""

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject)

    # search for a file with the correct name formatting
    file_pattern = f'{subject}_FinalGraph.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = np.loadtxt(run_file)

    return df

#not sure what this is
def get_objects(DFRAME):
    """ INPUT DataFrame 
       OUTPUT object sequence numbers for successor representation """
    data = drop_nan(DFRAME)
    obj_sequence = data["objnum"]
    return obj_sequence