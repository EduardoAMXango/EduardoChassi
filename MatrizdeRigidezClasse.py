import numpy as np
import matplotlib.pyplot as plt
import math

class Trelica():

    def __init__(self,area,elasticidade=1,comprimento=1,angulo=0):
        self.area= area
        self.elasticidade = elasticidade
        self.comprimento = comprimento
        self.angulo =  math.radians(angulo)
        
    def matriz_rigidez_nodal(self):
        c=(self.elasticidade*self.area)/(self.comprimento)
        K=np.array(([c ,-c],[-c, c]))
        print(K)   
        
    def matriz_rotacao(self):
        cos_teta=math.cos(self.angulo)
        sen_teta=math.sin(self.angulo)
        T = np.array([[cos_teta, sen_teta, 0, 0],
                    [-sen_teta, cos_teta, 0, 0],
                    [0, 0, cos_teta, sen_teta],
                    [0, 0, -sen_teta, cos_teta]])
        return T
    
    def matriz_rigidez_global(self): # Matriz de rigidez global do elemento de barra 2D
        cos_teta=math.cos(self.angulo)
        sen_teta=math.sin(self.angulo)
        c=(self.elasticidade*self.area)/(self.comprimento)
        K_e = c*np.array([[cos_teta**2, sen_teta*cos_teta, -cos_teta**2, -sen_teta * cos_teta],
                                   [0, cos_teta**2, -sen_teta*cos_teta, -sen_teta**2],
                                   [0, 0, cos_teta**2, sen_teta * cos_teta],
                                   [0, 0, 0, sen_teta**2]])
        K_e[np.tril_indices(4, k=-1)] = K_e.T[np.tril_indices(4, k=-1)]        # Preenchendo a parte inferior esquerda da matriz (simétrica)
        return K_e
  
viga = Trelica(8000,210,2000,60)
# Plotar a matriz de rotação
T = viga.matriz_rotacao()
plt.imshow(T, cmap='viridis', interpolation='nearest')
plt.title('Matriz de Rotação')
plt.colorbar()
plt.show()

# Plotar a matriz de rigidez global
K_global = viga.matriz_rigidez_global()
plt.imshow(K_global, cmap='viridis', interpolation='nearest')
plt.title('Matriz de Rigidez Global')
plt.colorbar()
plt.show()