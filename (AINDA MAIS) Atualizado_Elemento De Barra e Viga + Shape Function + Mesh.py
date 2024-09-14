import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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

    def calcular_comprimento(self, element):
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
        

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
        self.k_e = (E * I / L_e ** 3) * np.array([
            [12, 6 * L_e, -12, 6 * L_e],
            [6 * L_e, 4 * L_e ** 2, -6 * L_e, 2 * L_e ** 2],
            [-12, -6 * L_e, 12, -6 * L_e],
            [6 * L_e, 2 * L_e ** 2, -6 * L_e, 4 * L_e ** 2]
        ])

        self.m_e = (rho * A * L_e / 420) * np.array([
            [156, 22 * L_e, 54, -13 * L_e],
            [22 * L_e, 4 * L_e ** 2, 13 * L_e, -3 * L_e ** 2],
            [54, 13 * L_e, 156, -22 * L_e],
            [-13 * L_e, -3 * L_e ** 2, -22 * L_e, 4 * L_e ** 2]
        ])

        return self.k_e, self.m_e

    def matrizes_globais(self):
        for element in self.elements:
            node1, node2 = element
            self.k_e, self.m_e = self.matrizes_elementares(element)

            dofs = [2 * node1, 2 * node1 + 1, 2 * node2, 2 * node2 + 1]
            for i in range(4):
                for j in range(4):
                    self.K_global[dofs[i], dofs[j]] += self.k_e[i, j]
                    self.M_global[dofs[i], dofs[j]] += self.m_e[i, j]

        # Ancoragem
        ancoragem_nos = [0, 2, 4, 5]
        for no in ancoragem_nos:
            for i in range(self.num_dofs_per_node):
                idx = no * self.num_dofs_per_node + i
                self.K_global[idx, idx] = 1e20                          #Matriz de Rigidez GLO
                self.M_global[idx, idx] = 1e20                          #Matriz de Massa GLO

        return self.K_global, self.M_global

    def shape_fun(self, F1,F2,F3):
        E = 210e9   # Pa                                                #Módulo de Elasticidade
        I = 1.6667e-5 # m^4                                             #Momento de Inercia  
        G = 81.2e9  # Pa                                                #Módulo de Elasticidade em Cisalhamento
        A= 0.0125                                                       #Área da seção Transversal
        J = 1e-6    # m^4 (momento polar de inércia)
        torsao, deformacao, flexao = [], [], []
        for element in self.elements:
            L_e = self.calcular_comprimento(element)
            # Equação de torção
            torsao_val = (F1 * L_e) / (G * J)
            torsao.append(torsao_val)
            # Equação diferencial para deformação axial
            deformacao_val = (F2 / (A * E)) * L_e
            deformacao.append(deformacao_val)
            # Equação de Euler-Bernoulli para flexão
            x = np.linspace(0, L_e, len(F3))
            v_val = np.zeros_like(F3)
            for i in range(1, len(F3)):
                v_val[i] = (F3[i] - F3[i-1]) / (E * I) * (L_e ** 3 / 12)
            flexao.append(v_val)

        return np.array(torsao), np.array(deformacao), np.array(flexao)
    


    def Mesh(self):

        filename = input("Insira o nome do arquivo: ") + ".geo"
        diretorio = input("Insira o diretorio onde o arquivo .geo deve ser salvo: ")

        if not os.path.exists(diretorio):
            os.makedirs(diretorio)

        filepath = os.path.join(diretorio, filename)

        with open(filepath, 'w') as geo_file:
            for i, (x, y, z) in enumerate(self.nodes):
                geo_file.write(f'Point({i + 1}) = {{{x}, {y}, {z}, 1.0}};\n')

            for i, (start, end) in enumerate(self.elements):
                geo_file.write(f'Line({i + 1}) = {{{start + 1}, {end + 1}}};\n')

            if len(self.elements) > 2:
                line_loop_indices = ', '.join(str(i + 1) for i in range(len(elements)))
                geo_file.write(f'Line Loop(1) = {{{line_loop_indices}}};\n')
                geo_file.write('Plane Surface(1) = {1};\n')

            geo_file.write('Mesh.Algorithm = 6;\n')
            geo_file.write('Mesh.ElementOrder = 1;\n')
            geo_file.write('Mesh.Format = 1;\n')

        print(f'O arquivo foi salvo em: {filepath}, basta abrir o GMSH, e abrir o arquivo')

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
F1 = np.array([1000, 2000, 3000, 4000, 5000])
F2 = np.array([1000, 2000, 3000, 4000, 5000])
T = np.array([1000, 2000, 3000, 4000, 5000])
estrutura = Estrutura(elements, nodes, 1500, 8.33e-6, 8.33e-6)
K_global, M_global = estrutura.matrizes_globais()

#Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()


torsao, deformacao_axial, flexao  = estrutura.shape_fun(F1,F2,T)
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
axs[2].set_ylabel('Flexão(m)')
axs[2].legend()

plt.tight_layout()
plt.show()


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


