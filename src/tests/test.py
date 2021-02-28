import unittest
import gvar as gv
from meson_mass import Model

class MesonMassTest(unittest.TestCase):

    def test_sum(self):
        x = [0.03,  0.1,  0.4]
        m_D = ['1866.82(55)', '1876.81(97)', '1918.6(2.3)']
        m_B = ['5277.7(1.3)', '5288.2(1.3)', '5331.2(2.9)']
        m_Ds = ['1966.81(11)', '1968.87(25)', '1978.8(2.0)']
        m_Bs = ['5366.74(22)', '5369.82(41)', '5383.9(2.5)']

        self.assertEqual(list(map(str, Model('D').predict(x))), m_D,
                         "m_D({}) should be {}".format(x, m_D))
        self.assertEqual(list(map(str, Model('Ds').predict(x))), m_Ds,
                         "m_Ds({}) should be {}".format(x, m_Ds))
        self.assertEqual(list(map(str, Model('B').predict(x))), m_B,
                         "m_B({}) should be {}".format(x, m_B))
        self.assertEqual(list(map(str, Model('Bs').predict(x))), m_Bs,
                         "m_Bs({}) should be {}".format(x, m_Bs))

if __name__ == "__main__":
    unittest.main()
