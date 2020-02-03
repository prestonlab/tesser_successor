#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This util function is used to read in the behavioral data associated with TesserScan.
It gets the associated data directories for a given participant
    get_subj_dir(data_dir, subject_num)

It reads in the files associated with the following tasks:
- structured learning task
    load_struct_all(data_dir, subject_num)
    load_struct_run(data_dir, subject_num, part_num, run_num)
    drop_struct_nan(struct_data)
    get_struct_objects(struct_data)

- inductive inference task
    load_induct_

- grouping task

"""
import numpy as np
import pandas as pd
import os
from glob import glob


def get_subj_dir(data_dir, subject_num):
    """ INPUT data_dir, subject_num
       OUTPUT subject directory for a participant  """
    # check that the base directory exists
    if not os.path.exists(data_dir):
        raise IOError(f'Directory does not exist: {data_dir}')

    # look for directories with the correct pattern
    dir_search = glob(os.path.join(data_dir, f'tesserScan_{subject_num}_*'))
    if len(dir_search) != 1:
        raise IOError(f'Problem finding subject directory for {subject_num}')

    return dir_search[0]


def load_struct_df_all(data_dir, subject_num):
    """ INPUT data_dir, subject_num
       OUTPUT dataframe of all the structured learning  """

    # list of all runs to load
    parts = (1, 2)
    part_runs = {1: range(1, 6), 2: range(1, 7)}

    # load individual runs
    df_list = []
    for part in parts:
        for run in part_runs[part]:
            run_df = load_struct_run(data_dir, subject_num, part, run)
            df_list.append(run_df)

    # concatenate into one data frame
    df = pd.concat(df_list, sort=False)
    return df


def load_struct_df_run(data_dir, subject_num, part_num, run_num):
    """ INPUT data_dir, subject_num, part_num, run_num
       OUTPUT dataframe of a particular structured learning run  """

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject_num)

    # search for a file with the correct name formatting
    file_pattern = f'tesserScan_{subject_num}_*_StructLearn_Part{part_num}_Run_{run_num}.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject_num}, part {part_num}, run {run_num}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = pd.read_csv(run_file, sep='\t', skipinitialspace=True)

    # add a field indicating the experiment part
    df['part'] = part_num

    # remove trials with no object (these are just filler trials)
    df = df.loc[np.logical_not(np.isnan(df.objnum)), :]

    # convert object labels to integer
    df = df.astype({'objnum': 'int'})
    return df


def drop_struct_df_nan(struct_dframe):
    """ INPUT struct_dframe
       OUTPUT struct_dframe without NaNs """
    struct_dframe.replace(["NaN"], np.nan, inplace=True)
    data = struct_dframe.dropna()
    data = data.reset_index(drop=True)  # Resets the index to start at 0
    return data


def get_struct_objects(struct_dframe):
    """ INPUT struct_dframe
       OUTPUT just objects of struct_dframe """
    data = drop_nan(struct_dframe)
    obj_sequence = data["objnum"]
    return obj_sequence


def load_induct_df_all(data_dir, subject_num):
    """ INPUT data_dir, subject_num
       OUTPUT inductive generalization dataframe """

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject_num)

    # search for a file with the correct name formatting
    file_pattern = f'tesserScan_{subject_num}_*_InductGen.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject_num}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = pd.read_csv(run_file, sep='\t', skipinitialspace=True)
    return df


def load_induct_array_all(induct_dframe):
    """ INPUT induct_dframe
       OUTPUT four variable sequences for generalized induction """
    data = induct_dframe
    cue_sequence = data["CueNum"].values
    opt1_sequence = data["Opt1Num"].values
    opt2_sequence = data["Opt2Num"].values
    response_sequence = data["Resp"].values
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence


def load_group(data_dir, subject_num):
    """ INPUT data_dir, subject_num
       OUTPUT matrix of grouping data  """

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject_num)

    # search for a file with the correct name formatting
    file_pattern = f'{subject_num}_FinalGraph.txt'
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject_num}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    mat = np.loadtxt(run_file)
    return mat.astype(int)
