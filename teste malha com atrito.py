import pandas as pd
import matplotlib.pyplot as plt
from pyproj import Transformer
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.interpolate import griddata
import sys
import os
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(linewidth=200, suppress=True)

class CoordenadasUTM:
    def __init__(self, file_path, num_voltas=1, margem_malha=10, half_width=6.5):
        self.file_path = file_path
        self.transformer = Transformer.from_crs("EPSG:4326", "EPSG:32723", always_xy=True)  # UTM Zona 23S
        self.num_voltas = num_voltas
        self.margem_malha = margem_malha
        self.half_width = half_width

    def load_data(self):
        df = pd.read_excel(self.file_path, header=None, names=['Latitude', 'Longitude', 'Elevação'])
        utm = df.apply(lambda row: self.transformer.transform(row['Longitude'], row['Latitude']) + (row['Elevação'],), axis=1)
        self.utm_x, self.utm_y, self.utm_z = zip(*utm)
        min_x, min_y, min_z = min(self.utm_x), min(self.utm_y), min(self.utm_z)
        self.utm_x, self.utm_y, self.utm_z = [np.array(v) - min_v for v, min_v in zip([self.utm_x, self.utm_y, self.utm_z], [min_x, min_y, min_z])]

    def calc_track_edges(self):
        dx, dy = np.diff(self.utm_x), np.diff(self.utm_y)
        length = np.hypot(dx, dy)
        perpendicular = np.vstack([-dy, dx]) / length
        self.utm_x_left = self.utm_x[:-1] + perpendicular[0] * self.half_width
        self.utm_y_left = self.utm_y[:-1] + perpendicular[1] * self.half_width
        self.utm_x_right = self.utm_x[:-1] - perpendicular[0] * self.half_width
        self.utm_y_right = self.utm_y[:-1] - perpendicular[1] * self.half_width
        self.utm_x_left = np.append(self.utm_x_left, self.utm_x[-1] + perpendicular[0][-1] * self.half_width)
        self.utm_y_left = np.append(self.utm_y_left, self.utm_y[-1] + perpendicular[1][-1] * self.half_width)
        self.utm_x_right = np.append(self.utm_x_right, self.utm_x[-1] - perpendicular[0][-1] * self.half_width)
        self.utm_y_right = np.append(self.utm_y_right, self.utm_y[-1] - perpendicular[1][-1] * self.half_width)

    def generate_mesh(self):
        all_x, all_y = np.concatenate([self.utm_x, self.utm_x_left, self.utm_x_right]), np.concatenate([self.utm_y, self.utm_y_left, self.utm_y_right])
        min_x, max_x = min(all_x) - self.margem_malha, max(all_x) + self.margem_malha
        min_y, max_y = min(all_y) - self.margem_malha, max(all_y) + self.margem_malha
        grid_x, grid_y = np.meshgrid(np.linspace(min_x, max_x, 100), np.linspace(min_y, max_y, 100))
        grid_z = griddata((self.utm_x, self.utm_y), self.utm_z, (grid_x, grid_y), method='cubic')
        return grid_x, grid_y, grid_z

    def generate_friction_matrix(self, grid_x, grid_y):
        friction_matrix = np.full(grid_x.shape, 10)
        for i in range(len(self.utm_x) - 1):
            left = np.array([self.utm_x_left[i], self.utm_y_left[i]])
            right = np.array([self.utm_x_right[i], self.utm_y_right[i]])
            for j in range(grid_x.shape[0]):
                for k in range(grid_x.shape[1]):
                    point = np.array([grid_x[j, k], grid_y[j, k]])
                    if np.cross(right - left, point - left) < 1e-3:  # Se dentro da pista
                        friction_matrix[j, k] = 0.5
        return friction_matrix

    def animate_object(self, obj):
        def update(num):
            index = num % len(self.utm_x)
            obj.set_data([self.utm_x[index]], [self.utm_y[index]])
            obj.set_3d_properties([self.utm_z[index]])
        return update

    def plot(self):
        self.load_data()
        self.calc_track_edges()
        grid_x, grid_y, grid_z = self.generate_mesh()
        friction_matrix = self.generate_friction_matrix(grid_x, grid_y)
        print(friction_matrix)
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.utm_x_left, self.utm_y_left, self.utm_z, label='Extremidade Esquerda', color='green')
        ax.plot(self.utm_x_right, self.utm_y_right, self.utm_z, label='Extremidade Direita', color='red')

        ax.plot_surface(grid_x, grid_y, grid_z, facecolors=plt.cm.coolwarm(friction_matrix / 10), alpha=0.3, rstride=1, cstride=1, shade=False)

        obj, = ax.plot([self.utm_x[0]], [self.utm_y[0]], [self.utm_z[0]], 'bo')
        ani = animation.FuncAnimation(fig, self.animate_object(obj), frames=len(self.utm_x) * self.num_voltas, interval=100, repeat=False)

        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Elevação (m)')
        ax.set_title('Pista 3D com Atrito')
        plt.legend()
        plt.show()

# Executando o programa
file_path = 'C:\\Users\\dudua\\OneDrive\\Documentos\\GitHub\\EduardoChassi\\coordenadas.xlsx'
utm_converter = CoordenadasUTM(file_path, num_voltas=100)
utm_converter.plot()
