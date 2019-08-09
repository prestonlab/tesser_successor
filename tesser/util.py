#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# functions to read behavioral data (structured learning, induction)
import numpy as np
import pandas as pd
import os

p = '/home/rodrigo/Dropbox/tesser_successor/Data/'


def read_files(SUBJECT=0, TYPE='', PART=[1, 2], RUN=list(range(1, 7)), PATH=p):
    '''
        simple for loop to obtain all data from .txt files in a given directory
        Inputs:
        SUBJECT: Integer representing subject # for tesser behavioral studies
        TYPE: String describing the type of experimental data;
        'structured' for structured learning or 'induction' for gen induction
        PART: list containing the part numbers of the structured learning exp.
        RUN: list containing the run numbers of the structured learning exp.
        PATH: string describing the path taken to access tesser data
    '''
    data = {}
    s = 'tesserScan_' + str(100 + SUBJECT)
    for r, d, f in os.walk(PATH):    # r=root, d=directories, f = files
        for file in f:
            if TYPE == 'structured':
                for part in PART:
                    for run in RUN:
                        try:
                            if file.startswith(s) and file.endswith(
                              'StructLearn_Part%s_Run_%s.txt' % (part, run)):
                                df = pd.read_csv(
                                    os.path.join(r, file), sep="\t")
                                data[s + '_Part%s_Run_%s.txt' %
                                     (part, run)] = df
                        except ValueError:
                            pass
            if TYPE == 'induction':
                if file.startswith(s) and file.endswith('InductGen.txt'):
                    df = pd.read_csv(os.path.join(r, file), sep="\t")
                    data[s + '_InductGen.txt'] = df
        key = []
        for i in sorted(data.keys()):
            key.append(i)
    return data, key


def drop_nan(DATA):
    '''Drops NaN values from DataFrame'''
    DATA.replace([" NaN"], np.nan, inplace=True)
    DATA = DATA.dropna()
    DATA = DATA.reset_index(drop=True)  # Resets the index to start at 0
    return DATA


def get_objects(DFRAME):
    '''iNPUT DataFrame 
       OUTPUT object sequence numbers for successor representation'''
    data = drop_nan(DFRAME)
    obj_sequence = data[' objnum']
    return obj_sequence


def get_induction_data(DFRAME):
    '''iNPUT DataFrame 
       OUTPUT four variable sequences for generalized induction'''
    data = drop_nan(DFRAME)
    cue_sequence = np.array(data[[' CueNum']])
    opt1_sequence = np.array(data[[' Opt1Num']])
    opt2_sequence = np.array(data[[' Opt2Num']])
    response_sequence = np.array(data[[' Resp']])
    return cue_sequence, opt1_sequence, opt2_sequence, response_sequence


def print_keys(KEYS):
    print('Subject %s has %s keys for this data set' % (int(KEYS[0][11:14]), len(KEYS)))
    for i in KEYS:
        print(i)
