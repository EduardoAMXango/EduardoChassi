import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Estrutura:
    def __init__(self, elements, nodes, m, rho, A, E, I, Id, Ip):
        self.elements = elements                                            #Matriz de elementos conectados
        self.num_elements = len(elements)                                   #Número de elementos
        self.nodes = nodes                                                  #Matriz de nós com suas posiççoes
        self.num_nodes = len(nodes)                                         #Número total de nós
        self.massa = m                                                      #Massa do carro (Kg)
        self.momento_inercia_direcao = Id                                   #Momento de inércia em relação à direção (kg.m^2)
        self.momento_inercia_plano = Ip                                     #Momento de inércia em relação ao plano (kg.m^2)
        self.rho = rho                                                      #Massa específica do material do elemento (kg/m^3)
        self.A = A                                                          #Área da seção transversal de um elemento (m^2)
        self.I = I                                                          #Momento de inércia de área (m^4)
        self.E = E                                                          #Módulo de Young (Pa)
        self.num_dofs_per_node = 6                                                          #Graus de liberdade por nó
        self.num_dofs = self.num_nodes * self.num_dofs_per_node                             #Total de Graus de liberdade (gdls)
        self.K_global_barra = np.zeros((self.num_elements+1, self.num_elements+1))            #Matriz de rigidez global para barra
        self.M_global_barra = np.zeros((self.num_elements+1, self.num_elements+1))            #Matriz de massa global para barra
        self.K_global_viga = np.zeros((2*self.num_elements+2, 2*self.num_elements+2))       #Matriz de rigidez global para viga
        self.M_global_viga = np.zeros((2*self.num_elements+2, 2*self.num_elements+2))       #Matriz de massa global para viga
        self.num_modes = 20                                                                 #Número de modos de vibração a serem retornados


    def node_loc_matrix(self, node_tags, node_coord):
        # Número de Nós
        num_nodes = len(node_tags)

        # Gerando uma matrix número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)

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

        # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
        CM = np.zeros((num_connect, 3), dtype=int)

        # Preenchendo a matriz:
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2

        print("\n Conexão   1º Nó   2º Nó")
        print(CM)

    def matrizes_globais_barra(self):
        #Calculando as matrizes de rigidez e massa de cada elemento
        for element in self.elements:
            node1, node2 = element
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])     #Cálculo da distância entre nós a partir da geometria analítica
            k_e = (self.E * self.A / L_e) * np.array([[1, -1],
                                                      [-1, 1]])
            
            m_e = (self.rho * self.A * L_e / 6) * np.array([[2, 1],
                                                            [1, 2]])
        
            #Montando as matrizes globais
            dofs = [node1, node1+1, node2, node2+1]
            for i in range(2):
                for j in range(2):
                    self.K_global_barra[dofs[i], dofs[j]] += k_e[i, j]
                    self.M_global_barra[dofs[i], dofs[j]] += m_e[i, j]
        
        #Aplicando engastes nas extremidades 
        self.K_global_barra[0, 0] = 10**10
        self.K_global_barra[1, 1] = 10**10                             
        self.K_global_barra[-1, -1] = 10**10
        self.K_global_barra[-2, -2] = 10**10                           
        
        return self.K_global_barra, self.M_global_barra


    def matrizes_globais_viga(self):
        #Calculando as matrizes de rigidez e massa de cada elemento
        for element in self.elements:
            node1, node2 = element
            x1, y1, z1 = self.nodes[node1]
            x2, y2, z2 = self.nodes[node2]
            L_e = np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])     #Cálculo da distância entre nós a partir da geometria analítica
            k_e = (2 * self.E * self.I / L_e ** 3) * np.array([[6, -3 * L_e, -6, -3 * L_e],
                                                                [3 * L_e, 2 * L_e ** 2, -3 * L_e, L_e ** 2],
                                                                [-6, -3 * L_e, 6, 3 * L_e],
                                                                [-3 * L_e, L_e ** 2, 3 * L_e, 2 * L_e ** 2]])
            
            m_e = (self.rho * self.A * L_e / 420) * np.array([[156, 22 * L_e, 54, -13 * L_e],
                                                            [22 * L_e, 4 * L_e**2, 13 * L_e, -3 * L_e**2],
                                                            [54, 13 * L_e, 156, -22 * L_e],
                                                            [-13 * L_e, -3 * L_e**2, -22 * L_e, 4 * L_e**2]])
        
            #Montando as matrizes globais
            dofs = [2*node1, 2*node1+1, 2*node2, 2*node2+1]
            for i in range(4):
                for j in range(4):
                    self.K_global_viga[dofs[i], dofs[j]] += k_e[i, j]
                    self.M_global_viga[dofs[i], dofs[j]] += m_e[i, j]
        
        #Aplicando engastes nas extremidades 
        self.K_global_viga[0, 0] = 10**10                             
        self.K_global_viga[-1, -1] = 10**10                           
        
        return self.K_global_viga, self.M_global_viga
    
    def modal_analysis_viga(self):
        #Análise modal por resolução do problema de autovalor e autovetor
        unsorted_eigenvalues, unsorted_eigenvectors = eigh(self.K_global_viga, self.M_global_viga)

        #Frequências naturais (raiz quadrada dos autovalores)
        unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2*np.pi)            #Divisão por 2*pi para converter para hertz


        #Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
        sorted_indices = np.argsort(unsorted_frequencies)                           #Ordena as frequências em ordem crescente
        top_indices = sorted_indices[:self.num_modes]                               #Seleciona os índices dos primeiros n modos

        eigenvalues = np.array(unsorted_eigenvalues)[top_indices]                   #Filtra os primeiros n autovalores
        eigenvectors = np.array(unsorted_eigenvectors)[:, top_indices]                 #Filtra os primeiros n autovetores
        frequencies = np.array(unsorted_frequencies)[top_indices]                   #Filtra as primeiras n frequências

        return eigenvalues, eigenvectors, frequencies
    
    def modal_analysis_barra(self):
        #Análise modal por resolução do problema de autovalor e autovetor
        unsorted_eigenvalues, unsorted_eigenvectors = eigh(self.K_global_barra, self.M_global_barra)

        #Frequências naturais (raiz quadrada dos autovalores)
        unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2*np.pi)            #Divisão por 2*pi para converter para hertz


        #Tratando os dados (tomando apenas as 20 primeiras frequências naturais)
        sorted_indices = np.argsort(unsorted_frequencies)                           #Ordena as frequências em ordem crescente
        top_indices = sorted_indices[:self.num_modes]                               #Seleciona os índices dos primeiros n modos

        eigenvalues = np.array(unsorted_eigenvalues)[top_indices]                   #Filtra os primeiros n autovalores
        eigenvectors = np.array(unsorted_eigenvectors)[:, top_indices]                 #Filtra os primeiros n autovetores
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
estrutura = Estrutura(elements, nodes, 1500, 7850, 0.0125, 210e9, 8.33e-6, 8.33e-6, 8.33e-6)
K_global_viga, M_global_viga = estrutura.matrizes_globais_viga()
K_global_barra, M_global_barra = estrutura.matrizes_globais_barra()

#Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

#Gerar autovalores, autovetores e frequências naturais
autovalores_viga, autovetores_viga, frequencias_viga = estrutura.modal_analysis_viga()
autovalores_barra, autovetores_barra, frequencias_barra = estrutura.modal_analysis_barra()

print("\n Matriz de Rigidez Global para a estrutura de vigas")
print(K_global_viga)

print("\n Matriz de Massa Global para a estrutura de vigas")
print(M_global_viga)

print("\n Matriz de Rigidez Global para a estrutura de barras")
print(K_global_barra)

print("\n Matriz de Massa Global para a estrutura de barras")
print(M_global_barra)

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

#Exibindo as frequências naturais e modos de vibração da estrutura
print("\n Frequências Naturais (ω) da estrutura montada por vigas:")
print(frequencias_viga)

print("\n Frequências Naturais (ω) da estrutura montada por barras:")
print(frequencias_barra)

#Plotagem dos modos de vibração para a estrutura de vigas
for mode_idx in range(len(autovalores_viga)):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Modo {mode_idx + 1} Viga')

    for i, (x, y, z) in enumerate(nodes):
        ax.scatter(x, y, z, color='b', s=100)
        ax.text(x, y, z, f'  {i}', color='black', fontsize=8)

    for node1, node2 in elements:
        x = [nodes[node1][0], nodes[node2][0]]
        y = [nodes[node1][1], nodes[node2][1]]
        z = [nodes[node1][2], nodes[node2][2]]
        ax.plot(x, y, z, 'b--')

    mode_shape = autovetores_viga[:, mode_idx]
    displacements = np.zeros((len(nodes), 3))

    for j, (x, y, z) in enumerate(nodes):
        if j > 0 and j < len(nodes) - 1:
            displacements[j, 0] = mode_shape[2 * j]
            displacements[j, 1] = mode_shape[2 * j + 1]
            displacements[j, 2] = 0

    deformed_nodes = np.array(nodes) + 0.1 * displacements

    for node1, node2 in elements:
        x = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
        y = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
        z = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
        ax.plot(x, y, z, 'r-')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.tight_layout()
    plt.show()

#Plotagem dos modos de vibração para a estrutura de barras
for mode_idx in range(len(autovalores_barra)):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Modo {mode_idx + 1} Barra')

    for i, (x, y, z) in enumerate(nodes):
        ax.scatter(x, y, z, color='b', s=100)
        ax.text(x, y, z, f'  {i}', color='black', fontsize=8)

    for node1, node2 in elements:
        x = [nodes[node1][0], nodes[node2][0]]
        y = [nodes[node1][1], nodes[node2][1]]
        z = [nodes[node1][2], nodes[node2][2]]
        ax.plot(x, y, z, 'b--')

    mode_shape = autovetores_barra[:, mode_idx]
    displacements = np.zeros((len(nodes), 3))

    for j, (x, y, z) in enumerate(nodes):
        if j > 0 and j < len(nodes) - 1:
            displacements[j, 0] = mode_shape[j]
            displacements[j, 1] = mode_shape[j + 1]
            displacements[j, 2] = 0

    deformed_nodes = np.array(nodes) + 0.1 * displacements

    for node1, node2 in elements:
        x = [deformed_nodes[node1][0], deformed_nodes[node2][0]]
        y = [deformed_nodes[node1][1], deformed_nodes[node2][1]]
        z = [deformed_nodes[node1][2], deformed_nodes[node2][2]]
        ax.plot(x, y, z, 'r-')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.tight_layout()
    plt.show()

#PATCH NOTES DA NOVA VERSÃO (A 2.0):
#   -No init agora tomamos diretamente os vetores nodes e elements (posição dos nós + elementos) ao invés de inicializar ele na classe
#    estrutura nas linhas depois dele (se vcs verem a versão antiga fica claro o que eu tô falando);
#   -Consertados os problemas nas matrizes de M e K globais para viga que tavam dando problema para a análise modal;
#   -Implementadas as matrizes globais de barra, bem como sua análise modal;
#   -Adicionados os plots dos modos de vibração (porém ainda não tenho estabelecidas as restrições de movimento corretas);
#   -Quando forem jogar a parte do GMSH dentro me avisem pra eu poder acompanhar pra caso precise ajustar ou a função de matriz de
#    localização ou de conectividade;
#   -Ainda tem que consertar a enumeração dos pontos mas isso é extremamente desimportante pra falar a vdd :^|