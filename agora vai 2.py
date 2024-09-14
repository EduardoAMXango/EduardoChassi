import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Propriedades físicas e geométricas
L = 1.0
h = 0.01
b = 0.05
ro = 7800
E = 2.1e11
I = b * (h ** 3) / 12

# Propriedades do MEF
nel = 10
nnos = nel + 1
gdl_no = 2
gdl_global = nnos * gdl_no
n_mode = 3  # Número de modos a serem visualizados

Le = L / nel
Ae = b * h
me = Ae * ro

# Conectividades
matgdl = np.zeros((nnos, 2), dtype=int)
for i in range(nnos):
    matgdl[i, :] = [i * 2, (i * 2) + 1]

conect = np.zeros((nel, 4), dtype=int)
for i in range(nel):
    conect[i, :] = [i * 2, (i * 2) + 1, (i * 2) + 2, (i * 2) + 3]

# Matrizes globais
M = np.zeros((gdl_global, gdl_global))
K = np.zeros((gdl_global, gdl_global))

for i in range(nel):
    # Matriz de massa elementar
    Me = (Le * me / 420) * np.array([
        [156, 22 * Le, 54, -13 * Le],
        [22 * Le, 4 * Le ** 2, 13 * Le, -3 * Le ** 2],
        [54, 13 * Le, 156, -22 * Le],
        [-13 * Le, -3 * Le ** 2, -22 * Le, 4 * Le ** 2]
    ])

    # Matriz de rigidez elementar
    Ke = (E * I / Le ** 3) * np.array([
        [12, 6 * Le, -12, 6 * Le],
        [6 * Le, 4 * Le ** 2, -6 * Le, 2 * Le ** 2],
        [-12, -6 * Le, 12, -6 * Le],
        [6 * Le, 2 * Le ** 2, -6 * Le, 4 * Le ** 2]
    ])

    # Matriz de transferência
    TRANSF = np.zeros((4, gdl_global))
    for j in range(4):
        TRANSF[j, conect[i, j]] = 1

    # Montagem das matrizes globais
    K += TRANSF.T @ Ke @ TRANSF
    M += TRANSF.T @ Me @ TRANSF


# Função de análise modal
def modal_analysis(K, M, num_modes=20):
    unsorted_eigenvalues, unsorted_eigenvectors = eigh(K, M)
    unsorted_frequencies = np.sqrt(unsorted_eigenvalues) / (2 * np.pi)

    sorted_indices = np.argsort(unsorted_frequencies)
    top_indices = sorted_indices[:num_modes]

    eigenvalues = unsorted_eigenvalues[top_indices]
    eigenvectors = unsorted_eigenvectors[:, top_indices]
    frequencies = unsorted_frequencies[top_indices]

    return eigenvalues, eigenvectors, frequencies


# Análise modal
autovalores, autovetores, frequencias = modal_analysis(K, M)


# Visualização em 3D
def plot_mode_shape(mode_shape, mode_number):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x = np.linspace(0, L, nnos)
    y = mode_shape[::2]  # Considerando apenas deslocamentos verticais

    ax.plot(x, y, zs=0, zdir='z', label=f'Modo {mode_number + 1}')
    ax.scatter(x, y, zs=0, zdir='z')
    ax.set_xlabel('Comprimento da viga (m)')
    ax.set_ylabel('Deslocamento (m)')
    ax.set_zlabel('Modo de vibração')
    ax.set_title(f'Modo de Vibração {mode_number + 1}')

    plt.show()


# Plotando os modos de vibração
for i in range(n_mode):
    plot_mode_shape(autovetores[:, i], i)

print("Matriz de massa global (M):")
print(M)
print("\nMatriz de rigidez global (K):")
print(K)
print("\nFrequências Naturais (ω):")
print(frequencias)
print("\nAutovetores (Modos de Vibração):")
print(autovetores)
