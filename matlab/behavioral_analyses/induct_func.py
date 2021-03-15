#!/usr/bin/env python
from tesser import util 
from tesser import network
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy import stats
import matplotlib.pyplot as plt

def induct_perf_prop(data_dir, subject):
    induct_full = util.load_induct_subject(data_dir, subject)
    full_len = len(induct_full)
    #primary
    induct_prim = induct_full[induct_full["QuestType"]=="Prim"]
    prim_len = len(induct_prim)
    #bound1
    induct_bound1 = induct_full[induct_full["QuestType"]=="Bound1"]
    bound1_len = len(induct_bound1)
    #bound2
    induct_bound2 = induct_full[induct_full["QuestType"]=="Bound2"]
    bound2_len = len(induct_bound2)

    #get overall mean
    overall_acc_total = len(induct_full[induct_full["Acc"]==1])
    overall_acc_prop = overall_acc_total/full_len
    overall_mis_total = len(induct_full[induct_full["Acc"]==0])
    overall_mis_prop = overall_mis_total/full_len

    overall_bias = overall_acc_prop - overall_mis_prop
    
    #get primary mean (acc + RT)
    prim_acc_total = len(induct_prim[induct_prim["Acc"]==1])
    prim_acc_prop = prim_acc_total/prim_len
    prim_mis_total = len(induct_prim[induct_prim["Acc"]==0])
    prim_mis_prop = prim_mis_total/prim_len

    prim_bias = prim_acc_prop - prim_mis_prop

    #get bound1 mean (acc + RT)
    bound1_acc_total = len(induct_bound1[induct_bound1["Acc"]==1])
    bound1_acc_prop = bound1_acc_total/bound1_len
    bound1_mis_total = len(induct_bound1[induct_bound1["Acc"]==0])
    bound1_mis_prop = bound1_mis_total/bound1_len

    bound1_bias = bound1_acc_prop - bound1_mis_prop  
    
    #get bound2 mean (acc + RT)
    bound2_acc_total = len(induct_bound2[induct_bound2["Acc"]==1])
    bound2_acc_prop = bound2_acc_total/bound2_len
    bound2_mis_total = len(induct_bound2[induct_bound2["Acc"]==0])
    bound2_mis_prop = bound2_mis_total/bound2_len

    bound2_bias = bound2_acc_prop - bound2_mis_prop
    
    return overall_bias, prim_bias, bound1_bias, bound2_bias