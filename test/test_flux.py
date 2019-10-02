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

    def test_checkJetArgs(self):

        EPS = 1.0e-8
        Y0 = np.array([0.05, 1.0e53, 0.1, 0.4, 4, 0, 0, 0, 1.0, 2.2, 0.1,
                       0.01, 0.99, 1.0e28])
        Z0 = {'z': 0.5}

        Y = Y0.copy()
        Z = {}
        for k in Z0:
            Z[k] = Z0[k]

        models = [-2, -1, 0, 1, 2, 4]
        s = 0

        for m in models:
            self.assertIsNone(flux.checkJetArgs(m, s, *Y0, **Z))

            # theta_obs
            Y[0] = 0.0
            self.assertIsNone(flux.checkJetArgs(m, s, *Y, **Z))
            Y[0] = -0.1
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[0] = 0.5*np.pi + EPS
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[0] = Y0[0]

            # E0
            Y[1] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[1] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[1] = Y0[1]

            # theta_C
            Y[2] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[2] = -0.1
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[2] = 0.5*np.pi + EPS
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[2] = Y0[2]

            # theta_W
            if m != -1:
                Y[3] = 0.0
                self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
                Y[3] = -0.1
                self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
                Y[3] = 0.5*np.pi + EPS
                self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
                Y[3] = Y0[3]

            # b
            if m == 4:
                Y[4] = 0.0
                self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
                Y[4] = -1.0
                self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
                Y[4] = Y0[4]

            # L0
            Y[5] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[5] = Y0[5]

            # t_s
            Y[7] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[7] = Y0[7]

            # n0
            Y[8] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[8] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[8] = Y0[8]

            # p
            Y[9] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = 2.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = 1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = Y0[9]

            s = 2
            Y[9] = 2.0
            self.assertIsNone(flux.checkJetArgs(m, s, *Y, **Z))
            Y[9] = 1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[9] = Y0[9]
            s = 0

            # eps_e
            Y[10] = 1.0
            self.assertIsNone(flux.checkJetArgs(m, s, *Y, **Z))
            Y[10] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[10] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[10] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[10] = Y0[10]

            # eps_B
            Y[11] = 1.0
            self.assertIsNone(flux.checkJetArgs(m, s, *Y, **Z))
            Y[11] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[11] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[11] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[11] = Y0[11]

            # xi_N
            Y[12] = 1.0
            self.assertIsNone(flux.checkJetArgs(m, s, *Y, **Z))
            Y[12] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[12] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[12] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[12] = Y0[12]

            # dL
            Y[13] = 0.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[13] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Y[13] = Y0[13]

            # z
            Z['z'] = -1.0
            self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)
            Z['z'] = Z0['z']

        # specific for cone
        m = -2
        Y[3] = 0.5*Y[2]
        self.assertRaises(ValueError, flux.checkJetArgs, m, s, *Y, **Z)

    def test_checkCocoonArgs(self):

        EPS = 1.0e-8
        Y0 = np.array([10.0, 1.0, 1.0e53, 5, 1.0e-5, 0, 0, 0, 1.0, 2.2, 0.1,
                       0.01, 0.99, 1.0e28])
        Z0 = {'z': 0.5}

        Y = Y0.copy()
        Z = {}
        for k in Z0:
            Z[k] = Z0[k]

        models = [3]
        s = 0

        for m in models:
            self.assertIsNone(flux.checkCocoonArgs(m, s, *Y0, **Z))

            # u_max
            Y[0] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[0] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[0] = Y0[0]

            # u_min
            Y[1] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[1] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[1] = Y0[1]

            # Ei
            Y[2] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[2] = -0.1
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[2] = Y0[2]

            # Mej_solar
            Y[4] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s,
                              *Y, **Z)
            Y[4] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s,
                              *Y, **Z)
            Y[4] = Y0[4]

            # L0
            Y[5] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[5] = Y0[5]

            # t_s
            Y[7] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[7] = Y0[7]

            # n0
            Y[8] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[8] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[8] = Y0[8]

            # p
            Y[9] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[9] = 2.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[9] = 1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[9] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[9] = Y0[9]

            # eps_e
            Y[10] = 1.0
            self.assertIsNone(flux.checkCocoonArgs(m, s, *Y, **Z))
            Y[10] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[10] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[10] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[10] = Y0[10]

            # eps_B
            Y[11] = 1.0
            self.assertIsNone(flux.checkCocoonArgs(m, s, *Y, **Z))
            Y[11] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[11] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[11] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[11] = Y0[11]

            # xi_N
            Y[12] = 1.0
            self.assertIsNone(flux.checkCocoonArgs(m, s, *Y, **Z))
            Y[12] = 1.0 + EPS
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[12] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[12] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[12] = Y0[12]

            # dL
            Y[13] = 0.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[13] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Y[13] = Y0[13]

            # z
            Z['z'] = -1.0
            self.assertRaises(ValueError, flux.checkCocoonArgs, m, s, *Y, **Z)
            Z['z'] = Z0['z']


if __name__ == "__main__":
    unittest.main()
