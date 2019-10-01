import unittest
import numpy as np
import afterglowpy.flux as flux


class TestFlux(unittest.TestCase):

    def compareArrayEqualSingle(a1, a2):
        return (a1 == a2).all()

    def compareArrayTuple(self, func, argin, out):
        res = func(*argin)
        self.assertEqual(len(res), len(out))

        for i, a in enumerate(res):
            self.assertTrue((a == out[i]).all())

    def test_checkTNu(self):
        a1_5 = np.arange(5)
        a1_4 = np.arange(4)
        a2_45 = np.random.rand(4, 5)
        a2_46 = np.random.rand(4, 6)
        b2_45 = np.random.rand(4, 5)
        s1_1 = np.array([3.0])
        s1_4 = np.empty(4)
        s1_4[:] = 3.0
        s2_45 = np.empty((4, 5))
        s2_45[:] = 3.0

        # Check singleton args
        self.compareArrayTuple(flux.checkTNu, (s1_1, s1_1), (s1_1, s1_1))
        self.compareArrayTuple(flux.checkTNu, (3.0, s1_1), (s1_1, s1_1))
        self.compareArrayTuple(flux.checkTNu, (s1_1, 3.0), (s1_1, s1_1))
        self.compareArrayTuple(flux.checkTNu, (3.0, 3.0), (s1_1, s1_1))

        # Check correct 1d array args
        self.compareArrayTuple(flux.checkTNu, (a1_4, a1_4), (a1_4, a1_4))
        self.compareArrayTuple(flux.checkTNu, (s1_1, a1_4), (s1_4, a1_4))
        self.compareArrayTuple(flux.checkTNu, (a1_4, s1_1), (a1_4, s1_4))
        self.compareArrayTuple(flux.checkTNu, (3.0, a1_4), (s1_4, a1_4))
        self.compareArrayTuple(flux.checkTNu, (a1_4, 3.0), (a1_4, s1_4))

        # Check correct 2d array args
        self.compareArrayTuple(flux.checkTNu, (a2_45, b2_45), (a2_45, b2_45))
        self.compareArrayTuple(flux.checkTNu, (s1_1, a2_45), (s2_45, a2_45))
        self.compareArrayTuple(flux.checkTNu, (a2_45, s1_1), (a2_45, s2_45))
        self.compareArrayTuple(flux.checkTNu, (3.0, a2_45), (s2_45, a2_45))
        self.compareArrayTuple(flux.checkTNu, (a2_45, 3.0), (a2_45, s2_45))

        # Check incorrect args
        self.assertRaises(ValueError, flux.checkTNu, a1_4, a1_5)
        self.assertRaises(ValueError, flux.checkTNu, a1_4, a2_45)
        self.assertRaises(ValueError, flux.checkTNu, a2_45, a2_46)

    def test_checkThetaPhiTNu(self):
        a1_4 = np.arange(4)
        a1_5 = np.arange(5)
        a1_6 = np.arange(6)
        a2_45 = np.random.rand(4, 5)
        a2_46 = np.random.rand(4, 6)
        s1_1 = np.array([3.0])
        s1_4 = np.empty(4)
        s1_4[:] = 3.0
        s1_5 = np.empty(5)
        s1_5[:] = 3.0
        s1_6 = np.empty(6)
        s1_6[:] = 3.0
        s1_7 = np.empty(7)
        s1_7[:] = 3.0
        s2_45 = np.empty((4, 5))
        s2_45[:] = 3.0

        # Check singleton args
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, s1_1, s1_1, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, s1_1, s1_1, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, 3.0, s1_1, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, s1_1, 3.0, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, s1_1, s1_1, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, s1_1, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, s1_1, 3.0, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, s1_1, s1_1, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, 3.0, 3.0, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, 3.0, s1_1, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, s1_1, 3.0, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, 3.0, s1_1),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, s1_1, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, s1_1, 3.0, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, 3.0, 3.0, 3.0),
                               (s1_1, s1_1, s1_1, s1_1))

        # Check correct 1d array args
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, a1_4, a1_4, a1_4),
                               (a1_4, a1_4, a1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, a1_4, a1_4, a1_4),
                               (s1_4, a1_4, a1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, s1_4, a1_4, a1_4),
                               (a1_4, s1_4, a1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, a1_4, s1_4, a1_4),
                               (a1_4, a1_4, s1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, a1_4, a1_4, s1_4),
                               (a1_4, a1_4, a1_4, s1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, a1_4, a1_4, a1_4),
                               (s1_4, a1_4, a1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, 3.0, a1_4, a1_4),
                               (a1_4, s1_4, a1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, a1_4, 3.0, a1_4),
                               (a1_4, a1_4, s1_4, a1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, a1_4, a1_4, 3.0),
                               (a1_4, a1_4, a1_4, s1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a1_4, 3.0, 3.0, 3.0),
                               (a1_4, s1_4, s1_4, s1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, a1_4, 3.0, 3.0),
                               (s1_4, a1_4, s1_4, s1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, a1_4, 3.0),
                               (s1_4, s1_4, a1_4, s1_4))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, 3.0, a1_4),
                               (s1_4, s1_4, s1_4, a1_4))

        # Check correct 2d array args
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a2_45, a2_45, a2_45, a2_45),
                               (a2_45, a2_45, a2_45, a2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (s1_1, a2_45, a2_45, a2_45),
                               (s2_45, a2_45, a2_45, a2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a2_45, s2_45, a2_45, a2_45),
                               (a2_45, s2_45, a2_45, a2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a2_45, a2_45, s2_45, a2_45),
                               (a2_45, a2_45, s2_45, a2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a2_45, a2_45, a2_45, s2_45),
                               (a2_45, a2_45, a2_45, s2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (3.0, 3.0, a2_45, a2_45),
                               (s2_45, s2_45, a2_45, a2_45))
        self.compareArrayTuple(flux.checkThetaPhiTNu,
                               (a2_45, a2_45, 3.0, 3.0),
                               (a2_45, a2_45, s2_45, s2_45))

        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_4, a1_4, a1_5)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_4, a1_5, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_5, a1_4, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_5, a1_4, a1_4, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_6, a1_4, a1_4, a1_5)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_6, a1_5, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_5, a1_6, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_5, a1_4, a1_4, a1_6)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          3.0, a1_4, a1_4, a1_5)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, 3.0, a1_5, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_4, a1_5, 3.0, a1_4)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a1_5, a1_4, a1_4, 3.0)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a2_45, a2_45, a2_45, a2_46)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a2_45, a2_45, a2_46, a2_45)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a2_45, a2_46, a2_45, a2_45)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          a2_46, a2_45, a2_45, a2_45)
        self.assertRaises(ValueError, flux.checkThetaPhiTNu,
                          3.0, a2_45, a2_45, a1_4)


if __name__ == "__main__":
    unittest.main()
