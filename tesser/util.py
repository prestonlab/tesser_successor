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
from tesser import network


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

def load_struct_run_info(event_dir, subject, run):
    """Load event info for one structured learning run."""

    # search for a file with the correct name formatting
    file_pattern = f'tesser_{subject}_run{run}_info.txt'
    file_search = glob(os.path.join(event_dir, file_pattern))
    if len(file_search) != 1:
        raise IOError(f'Problem finding data for {subject}, run {run}.')
    run_file = file_search[0]

    # read log, fixing problem with spaces in column names
    run_info = pd.read_csv(run_file, sep='\,', skipinitialspace=True, header=None)
    
    #name the columns:
    run_info = run_info.rename({0:'trial_num', 1:'onset', 2:'TR', 3:'seq_type', 4:'item', 5:'dur'}, axis=1)
        
    return run_info

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
    nodes = network.temp_node_info()
    df['community'] = nodes.loc[df['objnum'], 'comm'].to_numpy()
    return df


def load_struct_bids(data_dir, subjects=None):
    """Load structure learning data for multiple subjects in BIDS format."""
    raw = load_struct(data_dir, subjects)
    orientation = {'cor': 'canonical', 'rot': 'rotated'}
    response = {'c': 'canonical', 'n': 'rotated'}
    df = pd.DataFrame(
        {
            'subject': raw['SubjNum'],
            'part': raw['part'],
            'run': raw['run'],
            'trial': raw['trial'],
            'trial_type': raw['seqtype'].astype('Int64'),
            'community': raw['community'],
            'object': raw['objnum'],
            'orientation': raw['orientnam'].map(orientation).astype('category'),
            'response': raw['resp'].map(response).astype('category'),
            'response_time': raw['rt'],
            'correct': raw['acc'].astype('Int64'),
        }
    )
    return df


def drop_struct_df_nan(struct_dframe):
    """Remove null trials (NaNs) at the beginning of structure task scans."""

    data = struct_dframe[pd.notnull(struct_dframe['obj'])]
    data = data.reset_index(drop=True)  # Resets the index to start at 0
    return data


def get_struct_objects(struct_dframe):
    """Series of just objects of from structured learning task."""

    obj_sequence = struct_dframe["objnum"]
    return obj_sequence


def object_count_run(this_run_info):
    """Getting count of each object 1-21 when inputting a structure learning dataframe."""
    this_run = this_run_info.reset_index(drop=True)
    all_count = []
    for item in range(1, 22):
        these_items = this_run[this_run["objnum"]==item]
        count = len(these_items)
        all_count.append(count)
    return all_count


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
    df.loc[:, 'response'] = df.Resp - 1
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
    nodes = network.temp_node_info()
    df['community'] = nodes.loc[df['cue'] + 1, 'comm'].to_numpy()
    return df


def load_induct_bids(data_dir, subjects=None):
    """Load induction data in BIDs format."""
    raw = load_induct(data_dir, subjects)
    trial_type = {'Prim': 'primary', 'Bound1': 'boundary1', 'Bound2': 'boundary2'}
    df = pd.DataFrame(
        {
            'subject': raw['SubjNum'],
            'trial': raw['TrialNum'],
            'trial_type': raw['QuestType'].map(trial_type).astype('category'),
            'environment': raw['Environment'],
            'community': raw['community'],
            'cue': raw['CueNum'],
            'opt1': raw['Opt1Num'],
            'opt2': raw['Opt2Num'],
            'response': raw['Resp'].astype('Int64'),
            'response_time': raw['RT'],
            'correct': raw['Acc'].astype('Int64'),
        }
    )
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

def score_induct(induct):
    """Score induction task data."""

    nodes = network.temp_node_info()
    for i, trial in induct.iterrows():
        correct_comm = nodes.loc[trial.CueNum, 'comm']
        answers = nodes.loc[nodes['comm'] == correct_comm]['node'].to_numpy()
        opts = [trial.Opt1Num, trial.Opt2Num]
        induct.loc[i, 'correct'] = np.nonzero(np.isin(opts, answers))[0][0]
    induct = induct.astype({'correct': int})
    return induct
