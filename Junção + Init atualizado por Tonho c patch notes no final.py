import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Estrutura:
    def __init__(self, n, m, rho, A, E, I, Id, Ip):
        self.num_nodes = n                                                  #Número total de nós
        self.massa = m                                                      #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                   #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                     #Momento de inércia em relação ao plano (kg.m^2)
        self.rho = rho                                                      #Massa específica do material do elemento (kg/m^3)
        self.A = A                                                          #Área da seção transversal de um elemento (m^2)
        self.I = I                                                          #Momento de inércia de área (m^4)
        self.E = E                                                          #Módulo de Young (Pa)
        self.num_dofs_per_node = 6                                          #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node             #Total de Graus de liberdade (gdls)
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))            #Matriz de rigidez global
        self.M_global = np.zeros((self.num_dofs, self.num_dofs))            #Matriz de massa global
        self.num_modes = 20                                                 #Número de modos de vibração a serem analizados


    def node_loc_matrix(self, node_tags, node_coord):
        # Número de Nós
        num_nodes = len(node_tags)

        # Gerando uma matrix número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((num_nodes, 4),
                                   dtype=float)  # Na primeira coluna, teremos os índices, nas seguintes, x y e z.

        # Preenchendo a matriz de zeros
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i]  # Número do nó na primeira coluna
            node_loc_matrix[i][1] = x  # Coordenada x na segunda coluna
            node_loc_matrix[i][2] = y  # Coordenada y na terceira coluna
            node_loc_matrix[i][3] = z  # Coordenada z na quarta coluna

        print("\n   Nó   x   y   z")
        print(node_loc_matrix)

    def connect_matrix(self):
        # Número de conexões
        num_connect = len(self.elements)

        # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
        CM = np.zeros((num_connect, 3), dtype=int)

        # Preenchendo a matriz:
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)

    def matriz_rigidez_global(self):
        for element in self.elements:
            node1, node2 = element
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])     #Cálculo da distância entre nós a partir da geometria analítica
            k_e = np.zeros((12, 12))
            coef = (2 * self.E * self.I / L_e ** 3)
            k_e[:4, :4] = coef * np.array([
                [6, -3 * L_e, -6, -3 * L_e],
                [3 * L_e, 2 * L_e ** 2, -3 * L_e, L_e ** 2],
                [-6, -3 * L_e, 6, 3 * L_e],
                [-3 * L_e, L_e ** 2, 3 * L_e, 2 * L_e ** 2]])
            k_e[6:10, 6:10] = k_e[:4, :4]
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs = dofs_node1 + dofs_node2
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += k_e[i, j]
            
            self.K_global[0,0] = 10**10                             #Engastes aplicados na matriz de rigidez global
            self.K_global[-1,-1] = 10**10                           #Engastes aplicados na matriz de rigidez global
        return self.K_global
    
    def matriz_massa_global(self):
        for node1, node2 in self.elements:
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])           #Cálculo da distância entre nós a partir da geometria analítica
            m_e = (self.rho * self.A * L_e / 420) * np.array([
                [156, 22 * L_e, 54, -13 * L_e],
                [22 * L_e, 4 * L_e**2, 13 * L_e, -3 * L_e**2],
                [54, 13 * L_e, 156, -22 * L_e],
                [-13 * L_e, -3 * L_e**2, -22 * L_e, 4 * L_e**2]
            ])

            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(4)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(4)]
            dofs = dofs_node1 + dofs_node2
            for i in range(4):
                for j in range(4):
                    self.M_global[dofs[i], dofs[j]] += m_e[i, j]
        return self.M_global
    
    def modal_analysis(self):
        #Análise modal por resolução do problema de autovalor e autovetor
        unsorted_eigenvalues, unsorted_eigenvectors = eigh(self.K_global, self.M_global)

        #Frequências naturais (raiz quadrada dos autovalores)
        unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2*np.pi)            #Divisão por 2*pi para converter para hertz


        #Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
        sorted_indices = np.argsort(unsorted_frequencies)                           #Ordena as frequências em ordem crescente
        top_indices = sorted_indices[:self.num_modes]                               #Seleciona os índices dos primeiros n modos

        eigenvalues = np.array(unsorted_eigenvalues)[top_indices]                   #Filtra os primeiros n autovalores
        eigenvectors = np.array(unsorted_eigenvectors)[top_indices]                 #Filtra os primeiros n autovetores
        frequencies = np.array(unsorted_frequencies)[top_indices]                   #Filtra as primeiras n frequências

        return eigenvalues, eigenvectors, frequencies

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
estrutura = Estrutura(len(nodes), 1500, 7850, 0.225, 210e9, 8.33e-6, 8.33e-6, 8.33e-6)
estrutura.nodes = nodes
estrutura.elements = elements
K_global = estrutura.matriz_rigidez_global()
M_global = estrutura.matriz_massa_global()

#Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

#Gerar autovalores, autovetores e frequências naturais
#autovalores, autovetores, frequencias = estrutura.modal_analysis()

print("\n Matriz de Rigidez Global")
print(K_global)

print("\n Matriz de Massa Global")
print(M_global)

# Convertendo a matriz de conectividade para um array numpy
connectivity_array = np.array(elements)

# Plotando o gráfico 3D da estrutura
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Adicionando os pontos
for i, (x, y, z) in enumerate(nodes):
    ax.scatter(x, y, z, color='b', s=100)
    ax.text(x, y, z, f'  {i}', color='black', fontsize=12)

# Adicionando as linhas de ligação entre os nós
for node1, node2 in elements:
    x = [nodes[node1][0], nodes[node2][0]]
    y = [nodes[node1][1], nodes[node2][1]]
    z = [nodes[node1][2], nodes[node2][2]]
    ax.plot(x, y, z, marker='o')

# Configurações adicionais
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Estrutura 3D')

plt.show()


#PATCH NOTES DA NOVA VERSÃO:
#   -Adicionados os outros parâmetros de geometria e material dos elementos (das barras do chassi) no init;
#   -Adicionado self.num_nodes no init (número de modos de vibração a serem retornados) e a inicialização da self.M_global;
#   -Generalizado a variável número de nós para que seja sempre um valor a partir da quantidade de nós no array nodes;
#   -Engastes aplicados na matriz de rigidez global;
#   -Coordenadas dos nós convertidas de milímetros para metros;
#   -Adicionada a parte de matriz de massa e ela está retornando a matriz de massa global;
#   -Adicionada a função de análise modal, PORÉM ela apresenta erro se utilizada (descomentada a linha 157 que pede os resultados);
#   -O erro provavelmente se dá devido à incompatibilidade das matrizes de massa e de rigidez global, ou seja, a forma de elas estarem sendo
#    montadas pode estar errado.