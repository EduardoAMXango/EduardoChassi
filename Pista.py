import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.interpolate import griddata

class CoordenadasUTM:
    def __init__(self, file_path):
        self.file_path = file_path
        self.df = None
        self.utm_coordinates = None
        self.utm_x = None
        self.utm_y = None
        self.utm_z = None
        self.utm_x_left = None
        self.utm_y_left = None
        self.utm_x_right = None
        self.utm_y_right = None
        self.transformer = self._get_transformer()
        self.num_voltas = 1  # Número de voltas
        self.margem_malha = 10  # Margem para estender a malha além da pista (em metros)

    # Método para carregar e processar os dados do Excel
    def load_and_process_data(self):
        self.df = pd.read_excel(self.file_path, header=None, names=['Latitude', 'Longitude', 'Elevação'])

    # Método privado para definir o transformador UTM
    def _get_transformer(self):
        return Transformer.from_crs("EPSG:4326", "EPSG:32723", always_xy=True)  # UTM Zona 23S (EPSG:32723)

    # Método para converter latitude, longitude e elevação para UTM
    def convert_to_utm(self):
        def convert_row(row):
            x, y = self.transformer.transform(row['Longitude'], row['Latitude'])
            return x, y, row['Elevação']
        
        self.utm_coordinates = self.df.apply(convert_row, axis=1)
    
    # Método para ajustar as coordenadas para que o menor ponto seja zero
    def adjust_coordinates(self):
        utm_x, utm_y, utm_z = zip(*self.utm_coordinates)
        
        min_x, min_y, min_z = min(utm_x), min(utm_y), min(utm_z)
        
        self.utm_x = [x - min_x for x in utm_x]
        self.utm_y = [y - min_y for y in utm_y]
        self.utm_z = [z - min_z for z in utm_z]

    # Método para calcular as extremidades da pista com base no centro
    def calculate_track_edges(self):
        half_width = 6.5
        utm_x_left = []
        utm_y_left = []
        utm_x_right = []
        utm_y_right = []

        for i in range(1, len(self.utm_x)):
            dx = self.utm_x[i] - self.utm_x[i - 1]
            dy = self.utm_y[i] - self.utm_y[i - 1]
            length = np.sqrt(dx**2 + dy**2)
            
            perpendicular_x = -dy / length
            perpendicular_y = dx / length

            utm_x_left.append(self.utm_x[i - 1] + perpendicular_x * half_width)
            utm_y_left.append(self.utm_y[i - 1] + perpendicular_y * half_width)
            utm_x_right.append(self.utm_x[i - 1] - perpendicular_x * half_width)
            utm_y_right.append(self.utm_y[i - 1] - perpendicular_y * half_width)

        utm_x_left.append(self.utm_x[-1] + perpendicular_x * half_width)
        utm_y_left.append(self.utm_y[-1] + perpendicular_y * half_width)
        utm_x_right.append(self.utm_x[-1] - perpendicular_x * half_width)
        utm_y_right.append(self.utm_y[-1] - perpendicular_x * half_width)

        self.utm_x_left, self.utm_y_left = utm_x_left, utm_y_left
        self.utm_x_right, self.utm_y_right = utm_x_right, utm_y_right

    # Função de animação para mover o objeto
    def animate_object(self, num_frames, obj):
        def update(num):
            pos_index = num % len(self.utm_x)  # A posição do objeto é calculada ao longo da pista
            obj.set_data([self.utm_x[pos_index]], [self.utm_y[pos_index]])
            obj.set_3d_properties([self.utm_z[pos_index]])

        return update

    # Método para gerar uma malha estendida
    def generate_mesh(self):
        # Adicionando margem nas coordenadas mínimas e máximas, considerando também as extremidades
        all_x = self.utm_x + self.utm_x_left + self.utm_x_right
        all_y = self.utm_y + self.utm_y_left + self.utm_y_right

        min_x, max_x = min(all_x) - self.margem_malha, max(all_x) + self.margem_malha
        min_y, max_y = min(all_y) - self.margem_malha, max(all_y) + self.margem_malha

        # Criando uma grade regular
        grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, 100), 
                                     np.linspace(min_y, max_y, 100))

        # Interpolando os valores z para os pontos da grade
        grid_z = griddata((self.utm_x, self.utm_y), self.utm_z, (grid_x, grid_y), method='cubic')

        return grid_x, grid_y, grid_z
    # Método para plotar os dados em 3D com movimento de um objeto
    def plot_utm_coordinates_with_movement(self, num_voltas):
        self.num_voltas = num_voltas
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plotando as extremidades da pista
        ax.plot(self.utm_x_left, self.utm_y_left, self.utm_z, label='Extremidade Esquerda', color='green')
        ax.plot(self.utm_x_right, self.utm_y_right, self.utm_z, label='Extremidade Direita', color='red')

        # Adicionando nós nas extremidades (esquerda e direita)
        ax.scatter(self.utm_x_left, self.utm_y_left, self.utm_z, color='green', marker='o', s=10, label='Nós Esquerda')
        ax.scatter(self.utm_x_right, self.utm_y_right, self.utm_z, color='red', marker='o', s=10, label='Nós Direita')

        # Criando um ponto que será o objeto a se mover na pista
        obj, = ax.plot([self.utm_x[0]], [self.utm_y[0]], [self.utm_z[0]], 'bo', label='Objeto')

        # Adicionando a malha estendida
        grid_x, grid_y, grid_z = self.generate_mesh()
        ax.plot_surface(grid_x, grid_y, grid_z, cmap='terrain', alpha=0.6)

        # Configurações do gráfico
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Objeto Percorrendo a Pista')

        # Animação do objeto
        num_frames = len(self.utm_x) * self.num_voltas
        ani = animation.FuncAnimation(fig, self.animate_object(num_frames, obj), frames=num_frames, interval=100, repeat=False)

        plt.legend()
        plt.show()

    # Método principal que organiza o fluxo de execução
    def process_and_plot(self, num_voltas=1):
        self.load_and_process_data()
        self.convert_to_utm()
        self.adjust_coordinates()
        self.calculate_track_edges()  # Calcula as extremidades da pista
        self.plot_utm_coordinates_with_movement(num_voltas)

# Executando o programa
file_path = 'C:\\Users\\dudua\\OneDrive\\Documentos\\GitHub\\EduardoChassi\\coordenadas.xlsx'
utm_converter = CoordenadasUTM(file_path)
utm_converter.process_and_plot(num_voltas=100)  # O objeto vai percorrer a pista 100 vezes
