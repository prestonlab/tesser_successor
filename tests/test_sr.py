import unittest
import numpy as np
import pandas as pd
from tesser import sr


class SRTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.n_state = 6
        self.env_step = [0, 1, 2, 1, 0, 1, 2, 3, 4, 5, 4, 3]
        self.SR_init = np.zeros((self.n_state, self.n_state))
        self.gamma = .9
        self.alpha = .5
        self.objnum = np.hstack((np.arange(1, 4), np.arange(4, 0, -1)))

    def test_sr_trials(self):
        SR = sr.run_experiment(self.env_step, self.gamma, self.alpha,
                               self.SR_init.copy(), self.n_state)
        expected = (
            np.array([[0.85125, 0.3375, 0., 0., 0., 0.],
                      [0.1125, 0.97625, 0.225, 0., 0., 0.],
                      [0., 0.1125, 0.75, 0., 0., 0.],
                      [0., 0., 0., 0.5, 0., 0.],
                      [0., 0., 0., 0.225, 0.75, 0.],
                      [0., 0., 0., 0., 0.225, 0.5]]))
        np.testing.assert_allclose(SR, expected)

    def test_sr_df(self):
        # construct a data frame with fake data
        df_list = []
        for part in 1, 2:
            for run in 1, 2:
                df_run = pd.DataFrame({'objnum': self.objnum, 'part': part,
                                       'run': run})
                df_list.append(df_run)
        df = pd.concat(df_list, 0, ignore_index=True)

        # run SR learning through all parts and runs
        sr_mats = sr.learn_sr(df, self.alpha, self.gamma)

        # check the state after each run
        expected = {(1, 1): np.array([[0.9 , 0.  , 0.  , 0.  ],
                                      [0.4 , 0.99, 0.  , 0.  ],
                                      [0.  , 0.4 , 0.99, 0.  ],
                                      [0.  , 0.  , 0.4 , 0.9 ]]),
                    (1, 2): np.array([[1.17, 0.45, 0.  , 0.  ],
                                      [0.53, 1.22, 0.04, 0.  ],
                                      [0.02, 0.54, 1.22, 0.04],
                                      [0.  , 0.02, 0.57, 1.17]]),
                    (2, 1): np.array([[1.26, 0.59, 0.02, 0.  ],
                                      [0.57, 1.29, 0.06, 0.  ],
                                      [0.03, 0.57, 1.28, 0.06],
                                      [0.  , 0.03, 0.63, 1.26]]),
                    (2, 2): np.array([[1.28, 0.64, 0.03, 0.  ],
                                      [0.58, 1.32, 0.07, 0.  ],
                                      [0.03, 0.59, 1.29, 0.07],
                                      [0.  , 0.03, 0.65, 1.28]])}
        for key, mat in expected.items():
            np.testing.assert_allclose(sr_mats[key], mat, atol=.01)


if __name__ == '__main__':
    unittest.main()
