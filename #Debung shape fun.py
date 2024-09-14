#Debung shape fun
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse.linalg import spsolve
from scipy.sparse import lil_matrix
import os

class Estrutura:
    def __init__(self, elements, nodes, m, Id, Ip):
        self.elements = elements
        self.num_elements = len(elements)
        self.nodes = nodes
        self.num_nodes = len(nodes)
        self.massa = m
        self.momento_inercia_direcao = Id
        self.momento_inercia_plano = Ip
        self.num_dofs_per_node = 6  
        self.num_dofs = self.num_nodes * self.num_dofs_per_node
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))
        self.num_modes = 20

    def calcular_comprimento(self, element):                            #Função auxiliar para cálculo de comprimento dos elementos
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
    

    def shape_fun(self, F_flexao,F_axial,F_torcao):
        E = 210e9   	#Modulo de Young (Pa)
        I = 1.6667e-5 	#Momento de inercia (m^4)
        G = 81.2e9  	#Modulo de Cisalhamento(Pa)
        A= 0.0125	    #Área da seção do elemento (m^2)	
        J = I/2     	#Momento polar de inércia (m^4) 
        torcao, deformacao, flexao = [], [], []
        for element in self.elements:
            L_e = self.calcular_comprimento(element)
            # Equação de torsão
            torcao_val = (F_torcao * L_e) / (G * J)         #Fonte[1]
            torcao.append(torcao_val)
            # Equação  para deformação axial
            deformacao_val = (F_axial* L_e / (A * E))       #Fonte[2]
            deformacao.append(deformacao_val)
            # Equação para flexão
            flexao_val = (F_flexao*L_e**3)/(3 * E * I)      #Fonte[3]
            flexao.append(flexao_val)
        # Plotando os resultados das deformações
        fig, axs = plt.subplots(3, 1, figsize=(12, 18))

        # Plot da Torção
        axs[0].plot(torcao, 'o-', label=[f'Força {F}N' for F in F_torcao])
        axs[0].set_title('Deformação por Torção de cada Elemento')
        axs[0].set_xlabel('Elemento')
        axs[0].set_ylabel('Torção (rad)')
        axs[0].legend()

        # Plot da Deformação Axial
        axs[1].plot(deformacao, 's-', label=[f'Força {F}N' for F in F_axial])
        axs[1].set_title('Deformação Axial de cada Elemento')
        axs[1].set_xlabel('Elemento')
        axs[1].set_ylabel('Deformação (m)')
        axs[1].legend()

        # Plot da Flexão
        axs[2].plot(flexao,'o-', label=[f'Força {F}N' for F in F_flexao])
        axs[2].set_title('Deformação por Flexão de cada Elemento')
        axs[2].set_xlabel('Posição ao longo do elemento')
        axs[2].set_ylabel('Flexão(m)')
        axs[2].legend()

        plt.tight_layout()
        plt.show()

        return np.array(torcao), np.array(deformacao), np.array(flexao)

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

#Criar a estrutura e montar as matrizes de rigidez e massa globais
#Dados: n = len(nodes), 
#       m = 1500 kg, 
#       rho = 7850 kg/m^3
#       A = 0.225 m^2
#       E = 210e9  # Módulo de elasticidade em Pa
#       I = 8.33e-6  # Momento de inércia em m^4
#       Ip = Id = 8.33e-6 kg.m^2
F_flexao = np.array([1000, 2000, 3000, 4000, 5000])
F_axial = np.array([1000, 2000, 3000, 4000, 5000])
F_torcao = np.array([1000, 2000, 3000, 4000, 5000])
estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)

estrutura.shape_fun(F_flexao,F_axial,F_torcao)


###REFERÊNCIAS
#[1]Resistência dos Materiais, Hibbeler 5º edição, capítulo 5.4,pág 167 - equação 5.15
#[2]Mecânica dos Materiais,Beer 5º edição, capítulo 2.8,pág 79 - equação 2.7
#[3]Mecânica dos Materiais,Beer 5º edição, capítulo 9.3,pág 555 - equação 9.11
