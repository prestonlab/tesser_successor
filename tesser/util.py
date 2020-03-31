#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This util function is used to read in the behavioral data associated
with TesserScan.
It gets the associated data directories for a given participant:
    get_subj_dir(data_dir, subject_num)

It reads in the files associated with the following tasks
- structured learning task:
    load_struct_all(data_dir, subject_num)
    load_struct_run(data_dir, subject_num, part_num, run_num)
    drop_struct_nan(struct_data)
    get_struct_objects(struct_data)

- inductive inference task:
    load_induct_df_all(data_dir, subject_num)
    load_induct_array_all(induct_dframe)

- grouping task:
    load_group(data_dir, subject_num)
"""
import numpy as np
import pandas as pd
import scipy.spatial.distance as sd
import os
from glob import glob


def subj_list():
    """Get IDs of included tesser participants."""
    participant_list = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
                        110, 111, 112, 113, 114, 115, 116, 117, 119, 120,
                        121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
                        131, 132, 133, 135, 136, 137, 138]
    return participant_list


def get_subj_dir(data_dir, subject_num):
    """Get the path to the data directory for a subject.

    data_dir : str
        Path to the base directory with subdirectories for each
        subject.

    subject_num : int
        Subject number without initials.
    """

    # check that the base directory exists
    if not os.path.exists(data_dir):
        raise IOError(f'Directory does not exist: {data_dir}')

    # look for directories with the correct pattern
    dir_search = glob(os.path.join(data_dir, f'tesserScan_{subject_num}_*'))
    if len(dir_search) != 1:
        raise IOError(f'Problem finding subject directory for {subject_num}')

    return dir_search[0]


def load_struct_run(data_dir, subject_num, part_num, run_num):
    """Load dataframe for one structured learning run."""

    # subject directory
    subj_dir = get_subj_dir(data_dir, subject_num)

    # search for a file with the correct name formatting
    file_pattern = (f'tesserScan_{subject_num}_*_StructLearn' 
                    f'_Part{part_num}_Run_{run_num}.txt')
    file_search = glob(os.path.join(subj_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject_num}, ' 
                      f'part {part_num}, run {run_num}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    df = pd.read_csv(run_file, sep='\t', skipinitialspace=True)

    # add a field indicating the experiment part
    df['part'] = part_num

    # remove the fixation trials in the runs (these are just filler trials)
    df = df[pd.notnull(df['objnum'])]

    # convert object labels to integer
    df = df.astype({'objnum': 'int'})
    return df


def load_struct_subject(data_dir, subject_num):
    """Load dataframe with structured learning task for one subject."""

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


def load_struct(data_dir, subjects=None):
    """Load structure learning data for all subjects."""

    if subjects is None:
        subjects = subj_list()

    df_all = []
    for subject in subjects:
        df_subj = load_struct_subject(data_dir, subject)
        df_all.append(df_subj)
    df = pd.concat(df_all, axis=0, ignore_index=True)
    return df


def drop_struct_df_nan(struct_dframe):
    """Remove null trials (NaNs) at the beginning of structure task scans."""

    data = struct_dframe[pd.notnull(struct_dframe['obj'])]
    data = data.reset_index(drop=True)  # Resets the index to start at 0
    return data


def get_struct_objects(struct_dframe):
    """Series of just objects of from structured learning task."""

#   data = drop_struct_df_nan(struct_dframe)
    obj_sequence = struct_dframe["objnum"]
    return obj_sequence


def load_induct_subject(data_dir, subject_num):
    """Load dataframe of inductive generalization task for one subject."""

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
    df.loc[:, 'cue'] = df.CueNum - 1
    df.loc[:, 'opt1'] = df.Opt1Num - 1
    df.loc[:, 'opt2'] = df.Opt2Num - 1
    df.loc[:, 'response'] = np.astype(df.Resp - 1, 'int')
    return df


def load_induct(data_dir, subjects=None):
    """Load induction data for all subjects."""

    if subjects is None:
        subjects = subj_list()

    df_all = []
    for subject in subjects:
        df_subj = load_induct_subject(data_dir, subject)
        df_all.append(df_subj)
    df = pd.concat(df_all, axis=0, ignore_index=True)
    return df


def load_induct_array_all(induct_dframe):
    """Unpack induction data into separate arrays.

    Returns arrays of cue_num, opt1_num, opt2_num, participant_resp.
    """

    data = induct_dframe
    cue_sequence = data["CueNum"].values
    opt1_sequence = data["Opt1Num"].values
    opt2_sequence = data["Opt2Num"].values
    response_sequence = data["Resp"].values
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence


def load_group_subject(data_dir, subject_num):
    """Load matrix of grouping data."""

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


# getting averages of matrix data, making sure that it is below diagonal etc.
def make_sym_matrix(asym_mat):
    """Calculate an average symmetric matrix from an asymmetric matrix."""

    v1 = sd.squareform(asym_mat, checks=False)
    v2 = sd.squareform(asym_mat.T, checks=False)
    vm = (v1 + v2) / 2
    sym_mat = sd.squareform(vm)
    return sym_mat
