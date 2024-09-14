import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import os

class Estrutura:
    def __init__(self, elements, nodes, m, Id, Ip):
        self.elements = elements                                             #Matriz de elementos conectados
        self.num_elements = len(elements)                                    #Número de elementos
        self.nodes = nodes                                                   #Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                          #Número total de nós
        self.massa = m                                                       #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                    #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                      #Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 2  # Correção: 2 graus de liberdade por nó para viga (deslocamento e rotação) #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node              #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))             #Matriz de massa global
        self.num_modes = 20                                                  #Número de modos de vibração a serem retornados

    def calcular_comprimento(self, element):                                 #Função auxiliar para cálculo de comprimento dos elementos
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def node_loc_matrix(self, node_tags, node_coord):
        num_nodes = len(node_tags)
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i] + 1
            node_loc_matrix[i][1] = x
            node_loc_matrix[i][2] = y
            node_loc_matrix[i][3] = z
        print("\n   Nó   x   y   z")
        print(node_loc_matrix)

    def connect_matrix(self):
        num_connect = len(self.elements)                                     # Número de conexões
        CM = np.zeros((num_connect, 3), dtype=int)
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2
        print("\n Conexão   1º Nó   2º Nó")
        print(CM)

    def matrizes_elementares(self, element):
        E = 210e9  # Pa
        I = 1.6667e-5  # m^4
        rho = 7850  # kg/m^3
        A = 0.0125  # m^2

        L_e = self.calcular_comprimento(element)
        k_e = (2 * E * I / L_e ** 3) * np.array([   [ 6       , -3 * L_e     , -6                  , -3 * L_e    , 0         , 0                    ],
                                                    [ 3 * L_e , 2 * L_e ** 2 , -3 * L_e            , L_e ** 2    , 0         , 0                    ],
                                                    [ 0       , 0            , (E * A / L_e)       , 0           , 0         , -(E * A / L_e)       ],
                                                    [ 0       , 0            , -6                  , -3 * L_e    , 6         , 3 * L_e              ],
                                                    [ 0       , 0            , -3 * L_e            , L_e ** 2    , 3 * L_e   , 2 * L_e ** 2         ],
                                                    [ 0       , 0            , -(E * A / L_e)      , 0           , 0         , (E * A / L_e)        ]])
        
        # Matriz de Massa Elementar
        m_e = (rho * A* L_e / 420) * np.array([     [156      , 22 * L_e     , 54                  , -13 * L_e   , 0         , 0                    ],
                                                    [22 * L_e , 4 * L_e**2   , 13 * L_e            , -3 * L_e**2 , 0         , 0                    ],
                                                    [ 0       , 0            , (rho * A * L_e / 3) , 0           , 0         , (rho * A * L_e / 6)  ],
                                                    [ 0       , 0            , 54                  , 13 * L_e    , 156       , -22 * L_e            ],
                                                    [ 0       , 0            , -13 * L_e           , -3 * L_e**2 , -22 * L_e , 4 * L_e**2           ],
                                                    [ 0       , 0            , (rho * A * L_e / 6) , 0           , 0         , (rho * A * L_e / 3)  ]])

        return k_e,m_e
    
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
    (0, 1),  # Elemento entre os nós 0 e 1 Nº0
    (1, 2),  # Elemento entre os nós 1 e 2 Nº1
    (4, 3),  # Elemento entre os nós 4 e 3 Nº2
    (3, 5),  # Elemento entre os nós 3 e 5 Nº3
    (1, 3)  # Elemento entre os nós 1 e 3  Nº4
]

#Criar a estrutura e montar as matrizes de rigidez e massa globais
#Dados: n = len(nodes), 
#       m = 1500 kg, 
#       rho = 7850 kg/m^3
#       A = 0.225 m^2
#       E = 210e9  # Módulo de elasticidade em Pa
#       I = 8.33e-6  # Momento de inércia em m^4
#       Ip = Id = 8.33e-6 kg.m^2

estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)

k_e,m_e=estrutura.matrizes_elementares(elements[3])

print("\n Matriz de Rigidez Global para a estrutura de vigas")
print(k_e)

print("\n Matriz de Massa Global para a estrutura de vigas")
print(m_e)