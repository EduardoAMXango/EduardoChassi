import numpy as np

class Estrutura:
    def __init__(self, E, I, nodes, elements):
        self.E = E
        self.I = I
        self.nodes = nodes
        self.elements = elements
        self.num_dofs_per_node = 6
        self.num_nodes = len(nodes)
        self.num_dofs = self.num_nodes * self.num_dofs_per_node
        self.K_global = np.zeros((self.num_dofs, self.num_dofs))

    def matriz_rigidez_elementar(self, L_e):
        # Matriz de rigidez elementar para um elemento com 6 DOFs por nó (total 12 DOFs)
        K_e = np.zeros((12, 12))
        coef = (2 * self.E * self.I / L_e**3)
        K_e[:4, :4] = coef * np.array([
            [6, -3 * L_e, -6, -3 * L_e],
            [3 * L_e, 2 * L_e**2, -3 * L_e, L_e**2],
            [-6, -3 * L_e, 6, 3 * L_e],
            [-3 * L_e, L_e**2, 3 * L_e, 2 * L_e**2]
        ])
        K_e[6:10, 6:10] = K_e[:4, :4]
        return K_e

    def calcular_comprimento(self, node1, node2):
        x1, y1, z1 = self.nodes[node1]
        x2, y2, z2 = self.nodes[node2]
        return np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)

    def mapear_dofs(self, element):
        node1, node2 = element
        dofs_node1 = [node1 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
        dofs_node2 = [node2 * self.num_dofs_per_node + i for i in range(self.num_dofs_per_node)]
        return dofs_node1 + dofs_node2

    def adicionar_rigidez(self, K_e, dofs):
        for i in range(len(dofs)):
            for j in range(len(dofs)):
                self.K_global[dofs[i], dofs[j]] += K_e[i, j]

    def montar_matriz_global(self):
        for element in self.elements:
            node1, node2 = element
            L_e = self.calcular_comprimento(node1, node2)
            K_e = self.matriz_rigidez_elementar(L_e)
            dofs = self.mapear_dofs(element)
            self.adicionar_rigidez(K_e, dofs)

        return self.K_global

# Exemplo de uso
E = 210e9  # Módulo de elasticidade em Pa
I = 8.33e-6  # Momento de inércia em m^4

L = 10
H = 1
W = 2

# Coordenadas dos nós (x, y, z)
nodes = [
    (0, 0, 0),
    (L, 0, 0),
    (L, H, 0),
    (0, H, 0),
    (0, 0, W),
    (L, 0, W),
    (L, H, W),
    (0, H, W)
]

# Conectividade dos elementos (índices dos nós)
elements = [
    (0, 1),  # Elemento entre os nós 0 e 1
    (1, 2),  # Elemento entre os nós 1 e 2
    (2, 3),  # Elemento entre os nós 2 e 3
    (4, 5),  # Elemento entre os nós 4 e 5
    (5, 6),  # Elemento entre os nós 5 e 6
    (6, 7)   # Elemento entre os nós 6 e 7
]

def node_loc_matrix(node_tags, node_coord):
    # Número de Nós
    num_nodes = len(node_tags)

    # Gerando uma matrix número de nos x 4: (Para cada linha, teremos: [índice do nó, x, y, z])
    node_loc_matrix = np.zeros((num_nodes, 4), dtype=float)  # Na primeira coluna, teremos os índices, nas seguintes, x y e z.

    # Preenchendo a matriz de zeros
    for i, (x, y, z) in enumerate(node_coord, start=0):
        node_loc_matrix[i][0] = node_tags[i]  # Número do nó na primeira coluna     
        node_loc_matrix[i][1] = x  # Coordenada x na segunda coluna
        node_loc_matrix[i][2] = y  # Coordenada y na terceira coluna
        node_loc_matrix[i][3] = z  # Coordenada z na quarta coluna

    print("   Nó   x   y   z")
    print(node_loc_matrix)

def connect_matrix(elements):
    # Número de conexões
    num_connect = len(elements)

    # Gerando a matriz de conectividade: (Para cada linha, teremos: [índice a coneção, 1º nó, 2º nó])
    CM = np.zeros((num_connect, 3), dtype=int)

    # Preenchendo a matriz:
    for i, (no1, no2) in enumerate(elements, start=0):
        CM[i][0] = i + 1
        CM[i][1] = no1
        CM[i][2] = no2

    print("Conexão   1º Nó   2º Nó")
    print(CM)

# Gerar as matrizes de localização dos nós e de conectividade
node_tags = list(range(len(nodes)))
node_loc_matrix(node_tags, nodes)
connect_matrix(elements)

# Criar a estrutura e montar a matriz de rigidez global
estrutura = Estrutura(E, I, nodes, elements)
K_global = estrutura.montar_matriz_global()

print("Matriz de Rigidez Global:")
print(K_global)
