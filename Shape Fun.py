import numpy as np
import matplotlib.pyplot as plt
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
        deform_axial = self.deformacao_axial(F_axial)
        flexao = self.flexao(F_flexao)
        torsao = self.torsao(F_torsao)
        self.print_results(deform_axial, flexao, torsao)
        return deform_axial, flexao, torsao

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

        # Condições de contorno
        for idx in [0, 1, -2, -1]:
            K_global[idx, :] = 0
            K_global[:, idx] = 0
            K_global[idx, idx] = 1e10 if idx % 2 else 1

        F_global = np.zeros(2 * num_nodes)
        for i in range(len(F_flexao)):
            F_global[2 * i] += F_flexao[i]

        v = spsolve(K_global.tocsr(), F_global)
        return v[::2]

    def torsao(self, F_torsao):
        return np.full(self.num_elements + 1, (F_torsao * self.L) / (self.G * self.J))

    def plot_results(self, deform_axial, flexao, torsao):
        x = np.linspace(0, self.L, self.num_elements + 1)

        plt.figure(figsize=(12, 8))

        # Plotar Deformação Axial
        plt.subplot(3, 1, 1)
        plt.plot(x, deform_axial, 'r-', marker='o')
        plt.title('Deformação Axial ao Longo da Estrutura')
        plt.xlabel('Comprimento (m)')
        plt.ylabel('Deformação Axial (m)')

        # Plotar Flexão
        plt.subplot(3, 1, 2)
        plt.plot(x, flexao, 'b-', marker='o')
        plt.title('Deformação por Flexão ao Longo da Estrutura')
        plt.xlabel('Comprimento (m)')
        plt.ylabel('Flexão (m)')

        # Plotar Torção
        plt.subplot(3, 1, 3)
        plt.plot(x, torsao, 'g-', marker='o')
        plt.title('Deformação por Torção ao Longo da Estrutura')
        plt.xlabel('Comprimento (m)')
        plt.ylabel('Torção (rad)')

        plt.tight_layout()
        plt.show()

    def print_results(self, deform_axial, flexao, torsao):
        print("Resultados da Deformação Axial:")
        print(deform_axial)
        print("\nResultados da Flexão:")
        print(flexao)
        print("\nResultados da Torção:")
        print(torsao)

# Exemplo de uso
A = 0.01    # m^2
E = 210e9   # Pa
I = 1.6667e-5 # m^4
G = 81.2e9  # Pa
K = 1000    # N/m
rho = 7850  # kg/m^3
g = 9.81    # m/s^2
J = 1e-6    # m^4 (momento polar de inércia)
L = 1.5     # m (comprimento total)
a = B = y = k = 1  # placeholder values
num_elements = 10

# Forças aplicadas
F_axial = 5000  # N
F_flexao = np.array([0, 2000, 3000, 4000, 5000])
F_torsao = 1000  # Nm

# Criar instância da classe e realizar a análise
estrutura = BodyAndFrameFEM(A, E, L, I, G, a, B, K, y, rho, k, g, J, num_elements)
deform_axial, flexao, torsao = estrutura.analyze(F_axial, F_flexao, F_torsao)

# Gerar os gráficos
estrutura.plot_results(deform_axial, flexao, torsao)