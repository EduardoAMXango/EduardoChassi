import numpy as np
import matplotlib.pyplot as plt

#Parâmetros físicos do ar:

p_atm= 101325                    #pressâo (N/m^2)
rho= 1.225                       #densidade(kg/m^3)
mi= 1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
ni= (mi/rho)                     #Viscosidade Cinematica (m²/s)
c_ar= 1000                       #Calor especifico (J/(kg*K))

#parâmetros do carro
length= 2                        #Comprimento do carro (m^2)
af= 1.5                          #Área Frontal do carro (m^2)
a_sup=2                          #Área de Superfície do carro (m^2)
cd_f = 0.2                       #Coeficiente de arrasto por atrito do carro
cd_p = 0.2                       #Coeficiente de arrasto por pressão do carro
ld= -0.3                         #Coeficiente de lift do carro

#parâmetros para a troca termica
deltaT=20                        #Variação de temperatura
a_res = 1                        #Área de resfriamento (m^2)
h= 200                           #Coeficiente de transferência de calor (W/m2*K)

#Equações Aerodinâmicas: 
def drag(drag_coeficient, pressao_dinamica,area):  # Função para o cálculo de arrasto
    return drag_coeficient * pressao_dinamica * area

def downforce(lift_coeficient,pressao_dinamica,area):  # Função para o cálculo de lift
    return lift_coeficient *pressao_dinamica*area

def numero_Reynolds(comprimento,velocidade,ni):     # Função para encontrar o numero de reynolds
    return (comprimento*velocidade)/ni
    
cd=cd_f+cd_p #Coeficiente de arrasto

# Velocidades de 0 a 80 m/s
velocidades = np.linspace(0.1,60,50)

for v in velocidades:
    # Cálculo da pressão dinâmica para cada velocidade
    pdinamica = 0.5 * rho  * (velocidades ** 2)

    # Cálculo do drag e lift para cada velocidade
    drags = drag(cd, pdinamica, af)
    downforces = downforce(ld, pdinamica,a_sup)

# Plotagem do gráfico
plt.figure(figsize=(10, 6))
plt.plot(velocidades, drags, label='Drag')
plt.plot(velocidades, downforces, label='Lift', linestyle='--')

# Configurações do gráfico
plt.title('Gráfico de Drag e Lift em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
plt.show()

#Arrefecimento
Qconvecção= h*a_res*deltaT #Calculo da quantidade de calor dissipado
print(f"{Qconvecção:0.2f} W")
