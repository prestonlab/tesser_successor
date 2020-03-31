import unittest
import numpy as np
from tesser import fit


class InductionCase(unittest.TestCase):

    def setUp(self) -> None:
        self.SR = np.array([[1, 2, 3],
                            [4, 5, 6],
                            [7, 8, 9]])

    def test_zero_mat(self):
        zero_sr = np.zeros((1, 3))
        tau = 1
        p = fit.probability_induction_choice(0, 1, 2, 0, zero_sr, tau)
        assert np.isnan(p)

    def test_choice_probability(self):
        cue = [0, 1, 2]
        opt1 = [1, 0, 1]
        opt2 = [2, 2, 0]
        response = [1, 0, 1]
        tau = 1
        p = [fit.probability_induction_choice(c, a, b, r, self.SR, tau) for
             c, a, b, r in zip(cue, opt1, opt2, response)]
        expected = [.6, .4, .4667]
        np.testing.assert_allclose(p, expected, atol=0.0001)

    def test_choice_by_tau(self):
        sr_row = np.array([[1, 2, 3]])
        p = [fit.probability_induction_choice(0, 1, 2, 0, sr_row, tau)
             for tau in [.5, 1, 2]]
        expected = [0.4494897427831781, 0.4, 0.3076923076923077]
        np.testing.assert_allclose(p, expected)


if __name__ == '__main__':
    unittest.main()
