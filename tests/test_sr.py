import pytest
import numpy as np
import pandas as pd
from tesser import sr


@pytest.fixture
def env_step():
    step = [0, 1, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3]
    return step


@pytest.fixture
def env_step_sr():
    expected = np.array(
        [[0.225, 0.85125, 0.1125, 0., 0., 0.],
         [0.25, 0.3375, 0.72625, 0., 0., 0.],
         [0., 0.25, 0.1125, 0.5, 0., 0.],
         [0., 0., 0., 0., 0.5, 0.],
         [0., 0., 0., 0.5, 0.225, 0.25],
         [0., 0., 0., 0., 0.5, 0.225]]
    )
    return expected


@pytest.fixture
def env_df():
    # construct a data frame with fake data
    objnum = np.hstack((np.arange(1, 4), np.arange(4, 0, -1)))
    df_list = []
    for part in 1, 2:
        for run in 1, 2:
            df_run = pd.DataFrame({'objnum': objnum, 'part': part,
                                   'run': run})
            df_list.append(df_run)
    df = pd.concat(df_list, 0, ignore_index=True)
    return df


@pytest.fixture
def env_df_sr():
    expected = {
        (1, 1): np.array(
            [[0., 0.9, 0., 0.],
             [0.9, 0.405, 0.09, 0.],
             [0., 0.9, 0.405, 0.09],
             [0., 0., 0.9, 0.405]]
        ),
        (1, 2): np.array(
            [[0.405, 1.17225, 0.0405, 0.],
             [1.09125, 0.5720625, 0.12735, 0.00405],
             [0.0405, 1.109475, 0.5356125, 0.12735],
             [0., 0.0405, 1.190475, 0.5315625]]
        ),
        (2, 1): np.array(
            [[0.5315625, 1.27465313, 0.0613575, 0.0018225],
             [1.15193813, 0.62924091, 0.14298694, 0.00659137],
             [0.0577125, 1.16332875, 0.57811978, 0.14116444],
             [0.0018225, 0.06217762, 1.28422125, 0.57152841]]
        ),
        (2, 2): np.array(
            [[0.57152841, 1.31062372, 0.07047987, 0.00314837],
             [1.17130423, 0.64842288, 0.1491612, 0.00783508],
             [0.06418313, 1.17832119, 0.59207482, 0.14601283],
             [0.00314837, 0.07115853, 1.31449232, 0.58423974]]
        ),
    }
    return expected


def test_sr_trials(env_step, env_step_sr):
    gamma = 0.9
    alpha = 0.5
    n_state = 6
    SR_init = np.zeros((n_state, n_state))
    SR = sr.run_experiment(env_step, gamma, alpha, SR_init.copy(), n_state)

    # TODO: check that this is actually the correct answer
    np.testing.assert_allclose(SR, env_step_sr)


def test_sr_df(env_df, env_df_sr):
    gamma = 0.9
    alpha = 0.5
    sr_mats = sr.learn_sr(env_df, alpha, gamma)
    for key, mat in env_df_sr.items():
        np.testing.assert_allclose(sr_mats[key], mat, atol=1e-6)
