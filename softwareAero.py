import numpy as np
import matplotlib.pyplot as plt

#Parâmetros físicos do ar:

p=1                             #pressâo (atm)
rho=1.184                       #densidade(kg/m^3)
mi=1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
ni=(mi/rho)                     #Viscosidade Cinematica (m²/s)

#parâmetros do carro
v=60                            #velocidade (m/s)
length=2                        #Comprimento do carro
a=0.6                           #Área Frontal do carro
drag_coeficient=0.3             #Coeficiente de arrasto do carro
lift_coeficient= 1              #Coeficiente de lift do carro

#Constantes para facilitar o calculo: 
denominador= (1/2 * rho*a*(v**2))

#Equações Aerodinâmicas: 

##num_Re= v*length/ni # Numero de Reynolds
drag= drag_coeficient*denominador # Equação do Arrasto
##lift2=lift_coeficient*denominador # Equação de lift

#Grafico de arrasto do corpo do carro:
velocidade = np.linspace(0.1,60, 5)  # Velocidade variando de 0 a 60 m/s
areas = np.linspace(0.2, 2, 5)  # Area variando de 0,2 a 2 m^2 
comprimento = 2 # Comprimento fixo em 2m

# Preparar a plotagem
plt.figure(figsize=(10, 6))

for area in areas:
    denominador2 = (1/2 * rho * area * (velocidade**2))
    num_Re2 = velocidade * comprimento * rho / mi  # Número de Reynolds
    dragf_coeficient = 0.031 / (num_Re2**(1/7)) # Calculo do coeficiente de arrasto causado pela fricção 
    dragt_coeficient=drag_coeficient + dragf_coeficient  #Somatorio dos coeficientes de arrasto
    drag = dragt_coeficient * denominador2  # Equação do Arrasto
    # Plotar o arrasto em função da velocidade para a área atual
    plt.plot(velocidade, drag, label=f'Área = {area:.2f} m²')

# Configurações do gráfico
plt.title('Gráfico de Drag em função da Velocidade para Diferentes Áreas')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Drag (N)')
plt.legend()
plt.grid(True)
plt.show()

