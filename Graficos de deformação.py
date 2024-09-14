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
        self.num_modes = 20                                                                 #Número de modos de vibração a serem retornados

    def calcular_comprimento(self, element):
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def elemento_barra_viga(self, F1, F2, T):
        torsao, deformacao, v = [], [], []
        for element in self.elements:
            L_e = self.calcular_comprimento(element)
            # Equação de torção
            torsao_val = (F1 * L_e) / (self.G * self.J)
            torsao.append(torsao_val)
            # Equação diferencial para deformação axial
            deformacao_val = (F2 / (self.A * self.E)) * L_e
            deformacao.append(deformacao_val)
            # Equação de Euler-Bernoulli para flexão
            x = np.linspace(0, L_e, len(T))
            v_val = np.zeros_like(T)
            for i in range(1, len(T)):
                v_val[i] = (T[i] - T[i-1]) / (self.E * self.I) * (L_e ** 3 / 12)
            v.append(v_val)

        return np.array(torsao), np.array(deformacao), np.array(v)

#Coordenadas dos nós (x, y, z)
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

F1 = np.array([1000, 2000, 3000, 4000, 5000])
F2 = np.array([1000, 2000, 3000, 4000, 5000])
T = np.array([1000, 2000, 3000, 4000, 5000])
estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)
torsao, deformacao_axial, flexao  = estrutura.elemento_barra_viga(F1,F2,T)
# Plotando os resultados das deformações
fig, axs = plt.subplots(3, 1, figsize=(12, 18))

# Plot da Torção
axs[0].plot(torsao, 'o-', label=[f'Elemento {x}' for x in range(5)])
axs[0].set_title('Deformação por Torção')
axs[0].set_xlabel('Elemento')
axs[0].set_ylabel('Torção (rad)')
axs[0].legend()

# Plot da Deformação Axial
axs[1].plot(deformacao_axial, 's-', label=[f'Elemento {x}' for x in range(5)])
axs[1].set_title('Deformação Axial')
axs[1].set_xlabel('Elemento')
axs[1].set_ylabel('Deformação (m)')
axs[1].legend()

# Plot da Flexão
for i, v_val in enumerate(flexao):
    axs[2].plot(v_val, label=f'Elemento {i}')
axs[2].set_title('Deformação por Flexão')
axs[2].set_xlabel('Posição ao longo do elemento')
axs[2].set_ylabel('Flexão')
axs[2].legend()

plt.tight_layout()
plt.show()