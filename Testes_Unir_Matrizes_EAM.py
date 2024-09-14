import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Estrutura:
    def __init__(self, A, E, I, G, K, rho, g, J, nodes, elements):
        self.A = A              # Área da seção transversal
        self.E = E              # Módulo de elasticidade
        self.I = I              # Momento de inércia
        self.G = G              # Módulo de elasticidade transversal
        self.K = K              # Constante de mola
        self.rho = rho          # Densidade
        self.g = g              # Aceleração devido à gravidade
        self.J = J              # Momento polar de inércia
        self.nodes = nodes        
        self.elements = elements
        self.num_dofs_per_node = 6
        self.num_nodes = len(nodes)
        self.num_dofs = self.num_nodes * self.num_dofs_per_node
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))

    def calcular_comprimento(self, element):
        node1, node2 = element
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def elemento_barra_viga(self, F1, F2, F3):
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
            x = np.linspace(0, L_e, len(F3))
            v_val = np.zeros_like(F3)
            for i in range(1, len(F3)):
                v_val[i] = (F3[i] - F3[i-1]) / (self.E * self.I) * (L_e ** 3 / 12)
            v.append(v_val)

        return np.array(torsao), np.array(deformacao), np.array(v)

    def connect_matrix(self):
        # Número de conexões
        num_connect = len(self.elements)
        # Gerando a matriz de conectividade
        CM = np.zeros((num_connect, 3), dtype=int)
        # Preenchendo a matriz
        for i, (no1, no2) in enumerate(self.elements, start=0):
            CM[i][0] = i + 1
            CM[i][1] = no1
            CM[i][2] = no2
        print("\n Conexão   1º Nó   2º Nó")
        print(CM)

    def node_loc_matrix(self, node_tags, node_coord):
        # Número de Nós
        num_nodes = len(node_tags)
        # Gerando uma matriz de localização dos nós
        node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)
        # Preenchendo a matriz de localização dos nós
        for i, (x, y, z) in enumerate(node_coord, start=0):
            node_loc_matrix[i][0] = node_tags[i]
            node_loc_matrix[i][1] = x
            node_loc_matrix[i][2] = y
            node_loc_matrix[i][3] = z
        print("\n   Nó   x   y   z")
        print(node_loc_matrix)

    def matriz_rigidez_global(self):
        for element in self.elements:
            node1, node2 = element
            L_e = self.calcular_comprimento(element)
            K_e = np.zeros((12, 12))
            coef = (2 * self.E * self.I / L_e**3)
            K_e[:4, :4] = coef * np.array([
                [6, -3 * L_e, -6, -3 * L_e],
                [3 * L_e, 2 * L_e**2, -3 * L_e, L_e**2],
                [-6, -3 * L_e, 6, 3 * L_e],
                [-3 * L_e, L_e**2, 3 * L_e, 2 * L_e**2]])
            K_e[6:10, 6:10] = K_e[:4, :4]
            dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
            dofs = dofs_node1 + dofs_node2
            for i in range(len(dofs)):
                for j in range(len(dofs)):
                    self.K_global[dofs[i], dofs[j]] += K_e[i, j]
        return self.K_global

# Exemplo de uso

A = 0.01    # m^2
E = 210e9   # Pa
I = 1.6667e-5 # m^4
G = 81.2e9  # Pa
K = 1000    # N/m
rho = 7850  # kg/m^3
g = 9.81    # m/s^2
J = 1e-6    # m^4 (momento polar de inércia)
F1 = np.array([1000, 2000, 3000, 4000, 5000])
F2 = np.array([1000, 2000, 3000, 4000, 5000])
F3 = np.array([1000, 2000, 3000, 4000, 5000])

# Coordenadas dos nós (x, y, z)
nodes = [
    (0, 0, 0),
    (0, 0.375, 0),
    (0, 0.700, 0),
    (1.500, 0.375, 0),
    (1.500, 0, 0),
    (1.500, 0.700, 0)
]

# Conectividade dos elementos (índices dos nós)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (4, 3),  # Elemento entre os nós 4 e 3
    (3, 5),  # Elemento entre os nós 3 e 5
    (1, 3)   # Elemento entre os nós 1 e 3
]

# Criar a estrutura e montar a matriz de rigidez global
estrutura = Estrutura(A, E, I, G, K, rho, g, J, nodes, elements)
K_global = estrutura.matriz_rigidez_global()

# Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
estrutura.node_loc_matrix(node_tags, nodes)
estrutura.connect_matrix()

print("\n Matriz de Rigidez Global")
print(K_global)

# Convertendo a matriz de conectividade para um array numpy
connectivity_array = np.array(elements)

# Calculando as deformações para torção, deformação axial e flexão
torsao, deformacao, v = estrutura.elemento_barra_viga(F1, F2, F3)

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

# Plotando os resultados das deformações
fig, axs = plt.subplots(3, 1, figsize=(12, 18))

# Plot da Torção
axs[0].plot(torsao, 'o-', label='Torção')
axs[0].set_title('Deformação por Torção')
axs[0].set_xlabel('Elemento')
axs[0].set_ylabel('Torção')
axs[0].legend()

# Plot da Deformação Axial
axs[1].plot(deformacao, 's-', label='Deformação Axial')
axs[1].set_title('Deformação Axial')
axs[1].set_xlabel('Elemento')
axs[1].set_ylabel('Deformação')
axs[1].legend()

# Plot da Flexão
for i, v_val in enumerate(v):
    axs[2].plot(v_val, label=f'Elemento {i}')
axs[2].set_title('Deformação por Flexão')
axs[2].set_xlabel('Posição ao longo do elemento')
axs[2].set_ylabel('Flexão')
axs[2].legend()

plt.tight_layout()
plt.show()