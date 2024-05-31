import numpy as np
import matplotlib.pyplot as plt
import math

#Declaração de Constantes:

p=1                             #pressâo (atm)
rho=1.184                       #densidade(kg/m^3)
mi=1.849*10**-5                 #Dynamic viscosity (kg/m*s)
ni=(mi/rho)                     #Kinematic viscosity (m²/s)

#Declaração de Variaveis 
v=50                            #velocidade (m/s)
length=2                        #Comprimento
width=1                         #Largura
a=width*length                  #Área

#Constantes para facilitar o calculo: 
denominador= (1/2 * rho*a*(v**2))

#Equações Aerodinâmicas: 

num_Re= v*length/ni # Numero de Reynolds
drag_coeficient = 1.328/(num_Re**0.5) # Equação do Coeficiente de Arrasto
drag= drag_coeficient*denominador # Equação do Arrasto

#Teste com Grafico:
velocidade = np.linspace(0.1,60, 5)  # Velocidade variando de 0 a 60 m/s
areas = np.linspace(0.2, 2, 5)  # Area variando de 0,2 a 2 m^2 
comprimento = 2 # Comprimento fixo em 2m

# Preparar a plotagem
plt.figure(figsize=(10, 6))

for area in areas:
    # Calcula o arrasto (Drag) para cada valor de velocidade
    denominador2 = (1/2 * rho * area * (velocidade**2))
    num_Re2 = velocidade * comprimento / ni  # Número de Reynolds
    drag_coeficient2 = 0.031 / (num_Re2**(1/7))  # Equação do Coeficiente de Arrasto
    drag2 = drag_coeficient2 * denominador2  # Equação do Arrasto

    # Plotar o arrasto em função da velocidade para a área atual
    plt.plot(velocidade, drag2, label=f'Área = {area:.2f} m²')
# Configurações do gráfico
plt.title('Gráfico de Drag em função da Velocidade para Diferentes Áreas')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Drag (N)')
plt.legend()
plt.grid(True)
plt.show()
print(drag)