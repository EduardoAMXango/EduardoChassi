import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
class Estrutura:
    def __init__(self, elements, nodes, m, Id, Ip):
        self.elements = elements                                            #Matriz de elementos conectados
        self.num_elements = len(elements)                                   #Número de elementos
        self.nodes = nodes                                                  #Matriz de nós com suas posiççoes
        self.num_nodes = len(nodes)                                         #Número total de nós
        self.massa = m                                                      #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                   #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                     #Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 6                                                          #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node                             #Total de Graus de liberdade (gdls)
        self.K_global_barra = np.zeros((self.num_elements+1, self.num_elements+1))            #Matriz de rigidez global para barra
        self.M_global_barra = np.zeros((self.num_elements+1, self.num_elements+1))            #Matriz de massa global para barra
        self.K_global_viga = np.zeros((2*self.num_elements+2, 2*self.num_elements+2))       #Matriz de rigidez global para viga
        self.M_global_viga = np.zeros((2*self.num_elements+2, 2*self.num_elements+2))       #Matriz de massa global para viga
        self.num_modes = 20   
                                                                      #Número de modos de vibração a serem retornados
    def Weak_Form(self):
        KF_total = 0
        KT_total = 0
        KF_elements = []
        KT_elements = []

        for elem in self.elements:
            n1, n2 = elem
            L = np.linalg.norm(self.nodes[n2] - self.nodes[n1])
            E = 210e9  # Módulo de elasticidade (Pa)
            G = 81.2e9  # Módulo de rigidez (Pa)
            A = 1e-4  # Área da seção transversal (m²)
            I = 1e-6  # Momento de inércia da seção transversal (m^4)
            J = 1e-6  # Momento polar de inércia (m^4)

            # Rigidez flexional
            KF = E * I / L
            # Rigidez torsional
            KT = G * J / L

            KF_total += KF
            KT_total += KT

            KF_elements.append(KF)
            KT_elements.append(KT)

        # Plotar os resultados
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))

        ax[0].bar(range(len(KF_elements)), KF_elements)
        ax[0].set_title('Rigidez Flexional por Elemento')
        ax[0].set_xlabel('Elemento')
        ax[0].set_ylabel('KF (N*m^2)')

        ax[1].bar(range(len(KT_elements)), KT_elements)
        ax[1].set_title('Rigidez Torsional por Elemento')
        ax[1].set_xlabel('Elemento')
        ax[1].set_ylabel('KT (N*m^2)')

        plt.tight_layout()
        plt.show()

        return KF_total, KT_total
nodes = [
    (0, 0, 0),
    (0, 0.375, 0),
    (0, 0.700, 0),
    (1.500, 0.375, 0),
    (1.500, 0, 0),
    (1.500, 0.700, 0)
]

#Conectividade dos elementos (índices dos nós)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (4, 3),  # Elemento entre os nós 4 e 3
    (3, 5),  # Elemento entre os nós 3 e 5
    (1, 3)  # Elemento entre os nós 1 e 3
]

#Criar a estrutura e montar as matrizes de rigidez e massa globais
#Dados: n = len(nodes), 
#       m = 1500 kg, 
#       rho = 7850 kg/m^3
#       A = 0.225 m^2
#       E = 210e9  # Módulo de elasticidade em Pa
#       I = 8.33e-6  # Momento de inércia em m^4
#       Ip = Id = 8.33e-6 kg.m^2
F1 = np.array([1000, 2000, 3000, 4000, 5000])
F2 = np.array([1000, 2000, 3000, 4000, 5000])
T = np.array([1000, 2000, 3000, 4000, 5000])
estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)
print(estrutura.Weak_Form())