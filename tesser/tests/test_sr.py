import unittest
import numpy as np
from tesser import sr


class SRTestCase(unittest.TestCase):

    def setUp(self) -> None:
        self.n_state = 6
        self.env_step = [1, 2, 3, 2, 1, 2, 3, 4, 5, 6, 5, 4]
        self.SR_init = np.zeros((self.n_state, self.n_state))
        self.gamma = .9
        self.alpha = .5

    def test_sr_learning(self):
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


if __name__ == '__main__':
    unittest.main()
