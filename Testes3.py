import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
class Estrutura:

    def __init__(self, elements, nodes, m, Id, Ip):
        self.elements = elements                                              # Matriz de elementos conectados
        self.num_elements = len(elements)                                     # Número de elementos
        self.nodes = nodes                                                    # Matriz de nós com suas posições
        self.num_nodes = len(nodes)                                           # Número total de nós
        self.massa = m                                                        # Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                     # Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                       # Momento de inércia em relação ao plano (kg.m^2)
        self.num_dofs_per_node = 6                                            # Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node               # Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))              # Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))              # Matriz de massa global
        self.num_modes = 12                                                   # Número de modos de vibração a serem retornados

    def node_loc_matrix(self, node_tags, node_coord):
        # Gerando uma matriz número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((self.num_nodes, 4), dtype=float)

        # Preenchendo a matriz de zeros
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i] + 1  # Número do nó na primeira coluna
            node_loc_matrix[i][1] = x  # Coordenada x na segunda coluna
            node_loc_matrix[i][2] = y  # Coordenada y na terceira coluna
            node_loc_matrix[i][3] = z  # Coordenada z na quarta coluna

        print("\n   Nó   x   y   z")
        print(node_loc_matrix)

    def connect_matrix(self):
        # Número de conexões
        num_connect = len(self.elements)

        # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a conexão, 1º nó, 2º nó])
        CM = np.zeros((num_connect, 3), dtype=int)

        # Preenchendo a matriz:
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)

    def calcular_comprimento(self, element):
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def element(self, elemento):
        # Variáveis e constantes físicas do modelo
        E = 210e9   # Pa
        I = 1.6667e-5 # m^4
        G = 81.2e9  # Pa
        A = 0.0125
        J = I / 2    # m^4 (momento polar de inércia)
        kappa = 0.5
        L_e = self.calcular_comprimento(elemento)
        Phi = (12 * E * I) / (kappa * G * A * L_e**2)
        Phi_barra = 1 / (1 + Phi)
        rho = 7850  # kg/m^3

        # Matriz de Rigidez Elementar
        k_e = np.array([
            [ E*A / L_e, 0,  0,                       0,                       0,                     0,            -E*A / L_e,            0,                       0,                       0,                       0,                     0],
            [ 0,        12 * Phi_barra * E*I / L_e**3,     0,                       0,                       6 * Phi_barra * E*I / L_e**2, 0,            0,               -12 * Phi_barra * E*I / L_e**3, 0,                       0,                       6 * Phi_barra * E*I / L_e**2, 0],
            [ 0,        0,  12 * Phi_barra * E*I / L_e**3,    0, 0,    -6 * Phi_barra * E*I / L_e**2, 0, 0,  -12 * Phi_barra * E*I / L_e**3, 0,   0, -6 * Phi_barra * E*I / L_e**2],
            [ 0,        0,                       0,                       G*J / L_e,                    0,                     0,            0,                       0,                       0,                       -G*J / L_e,                 0,                     0],
            [ 0,        6 * Phi_barra * E*I / L_e**2,     0,                       0,                       (4 + Phi) * Phi_barra * E*I / L_e, 0,            0,                       -6 * Phi_barra * E*I / L_e**2, 0,                       0,                       (2 - Phi) * Phi_barra * E*I / L_e, 0],
            [ 0,        0,                       -6 * Phi_barra * E*I / L_e**2,    0,                       0,                     (4 + Phi) * Phi_barra * E*I / L_e, 0,                       0,                       6 * Phi_barra * E*I / L_e**2, 0,                       0,                     (2 - Phi) * Phi_barra * E*I / L_e],
            [-E*A / L_e, 0,                       0,                       0,                       0,                     0,            E*A / L_e,                0,                       0,                       0,                       0,                     0],
            [ 0,        -12 * Phi_barra * E*I / L_e**3,  0,                       0,                       -6 * Phi_barra * E*I / L_e**2, 0,            0,                       12 * Phi_barra * E*I / L_e**3, 0,                       0,                       -6 * Phi_barra * E*I / L_e**2, 0],
            [ 0,        0,                       -12 * Phi_barra * E*I / L_e**3,  0,                       0,                     6 * Phi_barra * E*I / L_e**2, 0,                       0,                       12 * Phi_barra * E*I / L_e**3, 0,                       0,                     6 * Phi_barra * E*I / L_e**2],
            [ 0,        0,                       0,                       -G*J / L_e,                 0,                     0,            0,                       0,                       0,                       G*J / L_e,                  0,                     0],
            [ 0,        6 * Phi_barra * E*I / L_e**2,     0,                       0,                       (2 - Phi) * Phi_barra * E*I / L_e, 0,            0,                       -6 * Phi_barra * E*I / L_e**2, 0,                       0,                       (4 + Phi) * Phi_barra * E*I / L_e, 0],
            [ 0,        0,                       -6*Phi_barra*E*I/L_e**2,    0,                       0,                     (2 - Phi) * Phi_barra * E*I / L_e, 0,                       0,                       6 * Phi_barra * E*I / L_e**2, 0,                       0,                     (4 + Phi) * Phi_barra * E*I / L_e]
        ])
        return k_e

    def matrizes_global(self):
        #Calculando as matrizes de rigidez e massa de cada elemento
        for element in self.elements:
            el = element
            node1, node2 = element
            k_e = self.element(el)

            #Montando as matrizes globais
            #       x1  y1    z1    rx1   ry1   rz1   x2  y2    z2    rx2   ry2   rz2        
            #Montando as matrizes globais
            dofs = [6*node1, 6*node1+1,6*node1+2,6*node1+3,6*node1+4,6*node1+5, 6*node2, 6*node2+1,6*node2+2,6*node2+3,6*node2+4,6*node2+5]
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += k_e[i, j]
        
        #Aplicando engastes nas extremidades FAZER AS CONDIÇÕES DE CONTORNO NUMA DEF À PARTE
        self.K_global[0, 0] = 10**7
        self.K_global[1, 1] = 10**7                            
        self.K_global[-1, -1] = 10**7
        self.K_global[-2, -2] = 10**7                           
        
        return self.K_global
    
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

estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)
K_global= estrutura.matrizes_global()

k_elemental=estrutura.element(elements[0])
print(K_global)

print(len(K_global))