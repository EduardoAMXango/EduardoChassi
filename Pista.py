import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from mpl_toolkits.mplot3d import Axes3D

class Pista:
    def __init__(self, coordenadas, elevacao, atrito, largura_pista):
        self.gps_coords = coordenadas
        self.elevacao = elevacao
        self.atrito = atrito
        self.largura = largura_pista

    # Função para normalizar um vetor
    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
            return v
        return v / norm

    # Função para converter latitude/longitude para coordenadas planas X, Y
    def latlon_to_xy(self, lat, lon, lat0, lon0):
        R = 6371000  # Raio da Terra em metros
        x = R * np.radians(lon - lon0) * np.cos(np.radians(lat0))
        y = R * np.radians(lat - lat0)
        return x, y

    # Suavização da pista com interpolação cúbica
    def suavizacao(self):    
        lat0, lon0 = self.gps_coords[0]

        # Converter as coordenadas GPS para X, Y
        coordenadas_centro = [self.latlon_to_xy(lat, lon, lat0, lon0) + (1,) for lat, lon in self.gps_coords]

        # Separar as coordenadas X, Y, Z
        X_centro = [p[0] for p in coordenadas_centro]
        Y_centro = [p[1] for p in coordenadas_centro]
        Z_centro = [e/2 for e in self.elevacao]

        # Interpolação cúbica
        cs_x = CubicSpline(np.arange(len(X_centro)), X_centro, bc_type='periodic')
        cs_y = CubicSpline(np.arange(len(Y_centro)), Y_centro, bc_type='periodic')
        cs_z = CubicSpline(np.arange(len(Z_centro)), Z_centro, bc_type='periodic')

        # Gerar mais pontos
        t = np.linspace(0, len(X_centro)-1, 200)
        X_suave = cs_x(t)
        Y_suave = cs_y(t)
        Z_suave = cs_z(t)
        return X_suave, Y_suave, Z_suave

    # Cálculo das bordas esquerda e direita da pista
    def calcular_bordas(self, X, Y, Z):
        esquerda = []
        direita = []
        
        for i in range(len(X) - 1):
            vetor_direcao = np.array([X[i+1] - X[i], Y[i+1] - Y[i]])
            vetor_normal = np.array([-vetor_direcao[1], vetor_direcao[0]])
            vetor_normal = self.normalize(vetor_normal)
            
            esquerda.append((X[i] + vetor_normal[0] * self.largura, Y[i] + vetor_normal[1] * self.largura, Z[i]))
            direita.append((X[i] - vetor_normal[0] * self.largura, Y[i] - vetor_normal[1] * self.largura, Z[i]))
        
        esquerda.append((X[-1], Y[-1], Z[-1]))
        direita.append((X[-1], Y[-1], Z[-1]))
        
        return esquerda, direita

    # Geração da malha entre as bordas esquerda e direita
    def gerar_malha(self):
        X_suave, Y_suave, Z_suave = self.suavizacao()
        esquerda, direita = self.calcular_bordas(X_suave, Y_suave, Z_suave)

        X_esq, Y_esq, Z_esq = zip(*esquerda)
        X_dir, Y_dir, Z_dir = zip(*direita)

        X_malha = np.array([X_esq, X_dir])
        Y_malha = np.array([Y_esq, Y_dir])
        Z_malha = np.array([Z_esq, Z_dir])
        return X_malha, Y_malha, Z_malha, X_suave, Y_suave, Z_suave

    # Função para exibir a simulação 3D da pista
    def plotar_pista(self):
        X_malha, Y_malha, Z_malha, X_suave, Y_suave, Z_suave = self.gerar_malha()

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plotar o centro da pista
        ax.plot(X_suave, Y_suave, Z_suave, color='blue', label="Centro da Pista")

        # Separar as bordas
        X_esq, Y_esq, Z_esq = X_malha[0], Y_malha[0], Z_malha[0]
        X_dir, Y_dir, Z_dir = X_malha[1], Y_malha[1], Z_malha[1]

        # Plotar as bordas esquerda e direita
        ax.plot(X_esq, Y_esq, Z_esq, color='green', label="Borda Esquerda")
        ax.plot(X_dir, Y_dir, Z_dir, color='red', label="Borda Direita")

        # Conectar as bordas com linhas
        for i in range(len(X_esq)):
            ax.plot([X_esq[i], X_dir[i]], [Y_esq[i], Y_dir[i]], [Z_esq[i], Z_dir[i]], color='gray')

        # Plotar a malha
        ax.plot_surface(X_malha, Y_malha, Z_malha, color='black', alpha=0.7)

        # Configurações adicionais
        ax.set_xlabel('X (metros)')
        ax.set_ylabel('Y (metros)')
        ax.set_zlabel('Z (metros)')
        ax.set_title(f"Pista com Atrito de {self.atrito}")
        ax.legend()
        plt.show()

    # Verificar o ponto mais próximo da malha e exibir o atrito
    def verificar_atrito(self,ponto):
        """
        Verifica se o atrito está corretamente indexado para um ponto dado.
        
        Parâmetros:
        ponto (tuple): Um ponto (x, y, z) para verificar.
        X_malha, Y_malha, Z_malha (numpy arrays): Malha de pontos X, Y, Z.
        atrito (float): O valor do coeficiente de atrito aplicado.

        Retorna:
        dist_min (float): A menor distância encontrada para o ponto mais próximo.
        ponto_mais_proximo (tuple): O ponto da malha mais próximo do ponto dado.
        """
        x, y, z = ponto
        X_malha,Y_malha,Z_malha,_,_,_=self.gerar_malha()
        # Calcular a distância de cada ponto da malha até o ponto fornecido
        distancias = np.sqrt((X_malha - x) ** 2 + (Y_malha - y) ** 2 + (Z_malha - z) ** 2)
        
        # Encontrar a menor distância e o índice do ponto mais próximo
        indice_mais_proximo = np.unravel_index(np.argmin(distancias), distancias.shape)
        dist_min = distancias[indice_mais_proximo]
        
        # Obter as coordenadas do ponto mais próximo na malha
        ponto_mais_proximo = (X_malha[indice_mais_proximo], Y_malha[indice_mais_proximo], Z_malha[indice_mais_proximo])
        
        print(f"Coeficiente de atrito aplicado: {self.atrito}")
        print(f"Ponto mais próximo encontrado: {ponto_mais_proximo}")
        print(f"Distância mínima do ponto fornecido à malha: {dist_min:.4f} metros")
        
        return dist_min, ponto_mais_proximo

# Definir os dados da pista
gps_coords = [
    (-22.738762, -47.533146), (-22.739971, -47.533735), (-22.740344, -47.533928),
    (-22.740598, -47.533945), (-22.740725, -47.533782), (-22.740737, -47.533451),
    (-22.739432, -47.532570), (-22.739353, -47.531387), (-22.739159, -47.531069),
    (-22.738715, -47.530897), (-22.738259, -47.531082), (-22.737450, -47.531959),
    (-22.737394, -47.532273), (-22.737490, -47.532471), (-22.737608, -47.532600),
    (-22.738504, -47.533038), (-22.738762, -47.533146)  # Fechando o loop
]

elevacao = [10.4, 11.8, 11.7, 10.8, 9.8, 8.3, 7.5, 1.9, 0.4, 0.1, 1.4, 5.0, 6.6, 7.5, 8.0, 10.0, 10.4]

# Criar uma instância da classe Pista e gerar o gráfico
pista = Pista(gps_coords, elevacao, atrito=1.45, largura_pista=3)
pista.plotar_pista()

# Verificar o atrito em um ponto específico
ponto_teste = (228.1386, -2.8528, 1)
dist_min, ponto_mais_proximo = pista.verificar_atrito(ponto_teste)
X_malha,Y_malha,Z_malha,_,_,_=pista.gerar_malha()
