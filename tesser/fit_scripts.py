import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from . import util
from . import cython_sr
from . import cython_fit

from . import sr
from . import tasks
from . import network
from . import cfit
from . import csr
from scipy.spatial import distance
from scipy import optimize
from scipy.stats import linregress



def get_fits(
            struct_all, induct_all, fixed, 
            var_names, var_bounds, split,
             n_states, verbose, model_type,
             model, split_list
):
    results= cython_fit.fit_induct_indiv(
                                        struct_all, induct_all, fixed, 
                                         var_names, var_bounds, split,
                                         n_states, verbose, model_type,
                                         model, split_list
                                        )
if _init_ == '_main_':
    if is_fixed:
