import numpy as np
from scipy.stats import t as t_dist
import matplotlib.pyplot as plt

class Ridge():
    def __init__(self, X, T, r_A, r_B, r_g, U=None):
        """
           ([np.array]) X : a list of T(total time) by N dimensional abundance
                            matrix
           ([np.array]) U : a list of T(total time) by P dimensional perturbation
                            matrix; optional
           ([np.array]) T : a list of T by 1 dimensional number array containing
                            the timepoints at which measurements were taken
           (int) r_A : regularization parameter for interaction coefficients
           (int) r_B : regularization parameter for perturbation coefficients
           (int) r_g : regularization parameter for growth term
        """

        self.X = X
        self.U = U
        self.T = T
        self.r_A = r_A
        self.r_B = r_B
        self.r_g = r_g
        self.regularizers = self.make_regularizer_list()
        self.XU = []

        try:
            if self.U is not None:
                assert len(X) == len(U)
                for i in range(len(X)):
                    if U[i].shape[0] != X[i].shape[0]:
                        raise Exception("Error with dimensions.")
                    if np.sum([np.sum(u) for u in U]) != 0:
                        xu = np.hstack((np.exp(X[i]), np.ones((X[i].shape[0], 1)),
                            U[i]))
                        self.XU.append(xu[0:-1])
            else:
                for i in range(len(X)):
                    x = np.hstack((np.exp(X[i]), np.ones((X[i].shape[0], 1))))
                    self.XU.append(x[0:-1])

        except AssertionError:
            print("The dimensions of X and Y are not equal")

        self.Y = self.construct_Y()

    def construct_Y(self):
        """
            the left hand side of the gLV equation
        """

        Y = []
        for i in range(len(self.X)):
            X_i = self.X[i]
            shape = X_i.shape
            Y_arr = np.zeros((shape[0]-1, shape[1]))
            T_i = self.T[i]
            for t in range(1, shape[0]):
                del_t = T_i[t] - T_i[t-1]
                Y_arr[t-1] = (X_i[t] - X_i[t-1]) / del_t
            Y.append(Y_arr)

        return Y

    def make_regularizer_list(self):
        """
            combine the regularizers for different variables together
        """

        x_dim = self.X[0].shape[1]

        regularizers = []
        if self.U is None:
            regularizers = [self.r_A for i in range(x_dim)]+ [self.r_g]
        else:
            u_dim = self.U[0].shape[1]
            regularizers = [self.r_A for i in range(x_dim)] + [self.r_g] + \
                [self.r_B for i in range(u_dim)]

        return regularizers

    def solve(self):
        """
            solves the Ridge Regression problem and return the coefficients
        """

        x_dim = self.X[0].shape[1]
        Z_mat = np.vstack(self.XU)
        Y_mat = np.vstack(self.Y)
        diag_reg = np.diag(np.array(self.regularizers))

        ZZ = np.matmul(Z_mat.T, Z_mat)
        term1 = np.linalg.inv(ZZ + diag_reg)
        term2 = np.matmul(Z_mat.T, Y_mat)

        theta = np.matmul(term1, term2)
        A = theta[0:x_dim, :]
        g = theta[x_dim: x_dim + 1, :]
        B = None
        if self.U is not None:
            B = theta[x_dim+1:, :]

        return A, g, B

    def map_t_to_p(self, t_values, deg_freedom):
        """
            generates p-values associated with the test statistics for Student-t
            distribution
        """

        dist = t_dist(deg_freedom)
        p_values = 2 * (1 - dist.cdf(np.abs(t_values)))

        return p_values

    def compute_hat_matrix(self):
        """
            computes the hat matrix used in computation of effective degrees
            of freedom
        """

        Z_mat = np.vstack(self.XU)
        Y_mat = np.vstack(self.Y)
        diag_reg = np.diag(np.array(self.regularizers))

        ZZ = np.matmul(Z_mat.T, Z_mat)
        term1 = np.linalg.inv(ZZ + diag_reg)
        hat_mat = np.matmul(Z_mat, np.matmul(term1, Z_mat.T))

        return hat_mat

    def significance_test(self):
        """
            performs a significance test proposed in Lule et al. 2011;
            spits out a matrix of p-values associated with beta coefficients
            of the regression model
        """

        def compute_V_mat(Z, Y, D):
            """matrix used in computation of variance of beta coefficients"""

            mat1 = np.linalg.inv(np.matmul(Z.T, Z) + D)
            mat2 = np.matmul(Z.T, Z)
            V = np.matmul(mat1, np.matmul(mat2, mat1))

            return V

        Z_mat = np.vstack(self.XU)
        Y_mat = np.vstack(self.Y)
        diag_reg = np.diag(np.array(self.regularizers))
        theta_A, theta_g, theta_B = self.solve()

        try:
            Theta = np.vstack((theta_A, theta_g, theta_B))
        except ValueError:
            print("No valid perturbation present. Concatenating growth and"\
                "interaction data")
            Theta = np.vstack((theta_A, theta_g))

        n, m = Z_mat.shape
        H = self.compute_hat_matrix()
        nu_sigma = n - np.trace(2 * H - np.matmul(H, H.T))
        V_mat = compute_V_mat(Z_mat, Y_mat, diag_reg)

        diff_Y = Y_mat - np.matmul(Z_mat, Theta)
        Sigma = np.matmul(diff_Y.T, diff_Y)
        sigma = np.diagonal(Sigma) / nu_sigma

        t_stats = np.zeros(Theta.shape)
        for i in range(Y_mat.shape[1]):
            s = sigma[i]
            se = np.sqrt(s * np.diag(V_mat))
            theta_i = Theta[:, i]
            t_stats[:, i] = theta_i / se

        nu_t_test = n - np.trace(H)
        p_vals = self.map_t_to_p(t_stats, nu_t_test)
        x_dim = self.X[0].shape[1]

        p_B = None
        if self.U is not None:
            p_B = p_vals[x_dim+1:, :]
        p_A = p_vals[0:x_dim, :]
        p_g = p_vals[x_dim: x_dim + 1, :]

        return p_A, p_g, p_B
