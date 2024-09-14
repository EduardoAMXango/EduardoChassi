import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class BodyAndFrameFEM:
    def __init__(self, A, E, L, I, G, a, B, K, y, rho, k, g, J, num_elements):
        self.A = A
        self.E = E
        self.L = L
        self.I = I
        self.G = G
        self.a = a
        self.B = B
        self.K = K
        self.y = y
        self.rho = rho
        self.k = k
        self.g = g
        self.J = J
        self.num_elements = num_elements

    def analyze(self, F_axial, F_flexao, F_torsao):
        return (
            self.deformacao_axial(F_axial),
            self.flexao(F_flexao),
            self.torsao(F_torsao)
        )

    def deformacao_axial(self, F_axial):
        return np.full(self.num_elements + 1, (F_axial / (self.A * self.E)) * (self.L / self.num_elements))

    def flexao(self, F_flexao):
        num_nodes = self.num_elements + 1
        L_e = self.L / self.num_elements
        K_global = lil_matrix((2 * num_nodes, 2 * num_nodes))

        for i in range(self.num_elements):
            K_e = (self.E * self.I / L_e**3) * np.array([
                [12, 6*L_e, -12, 6*L_e],
                [6*L_e, 4*L_e**2, -6*L_e, 2*L_e**2],
                [-12, -6*L_e, 12, -6*L_e],
                [6*L_e, 2*L_e**2, -6*L_e, 4*L_e**2]
            ])
            dofs = [2*i, 2*i+1, 2*i+2, 2*i+3]
            for ii in range(4):
                for jj in range(4):
                    K_global[dofs[ii], dofs[jj]] += K_e[ii, jj]

        # condições de contorno
        for idx in [0, 1, -2, -1]:
            K_global[idx, :] = 0
            K_global[:, idx] = 0
            K_global[idx, idx] = 1e10 if idx % 2 else 1

        F_global = np.zeros(2 * num_nodes)
        for i in range(len(F_flexao)):
            F_global[2 * i] += F[i]

        v = spsolve(K_global.tocsr(), F_global)
        return v[::2]

    def torsao(self, F_torsao):
        return np.full(self.num_elements + 1, (F_torsao * self.L) / (self.G * self.J))





