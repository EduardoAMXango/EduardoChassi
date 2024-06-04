import numpy as np
import matplotlib.pyplot as plt

#Parâmetros físicos do ar:

p_atm= 101325                    #pressâo (N/m^2)
rho= 1.184                       #densidade(kg/m^3)
mi= 1.849*10**-5                 #Viscosidade dinâmica (kg/m*s)
ni= (mi/rho)                     #Viscosidade Cinematica (m²/s)
c_ar= 1                          #Calor especifico (kJ/(kg*K))
temp_ar= 25                      #Temperatura inicial do ar(°C)
temp_f=60

#parâmetros do carro
#v=60                            #velocidade (m/s)
length= 2                        #Comprimento do carro
af= 0.5                          #Área Frontal do carro
a_sup=2                          #Área de Superfície do carro
a1=0.25                          #Área de entrada do ar embaixo do carro
a2=0.20                          #Área embaixo do carro
a_cooling= 0.05                  #Área da entrada de ar para arrefecimento
cd = 0.75                         #Coeficiente de arrasto do carro
ld= -3                         #Coeficiente de lift do carro

#Equações Aerodinâmicas: 
def drag(drag_coeficient, pressao_dinamica,area):  # Função para o cálculo de arrasto
    return drag_coeficient * pressao_dinamica * area

def drag_friction(velocidade,rho,mi,comprimento,area):  # Função para o cálculo de arrasto causado pelo atrito
    numrey=velocidade*length/ni
    cd_friction= 0.074/(numrey**(1/5))
    return cd_friction*velocidade**2*area*0.5*rho

def lift(lift_coeficient,pressao_dinamica,area):  # Função para o cálculo de lift
    return lift_coeficient *pressao_dinamica*area

def lift_coeficient(lift,pressao_dinamica,area):  # Função para o cálculo do coeficiente de lift
    return lift/(pressao_dinamica*area)

# Velocidades de 0 a 60 m/s
velocidades = np.linspace(0.1,60,10)

# Cálculo da pressão dinâmica para cada velocidade
pdinamica = 0.5 * rho  * (velocidades ** 2)

#Bernoulli
v2=velocidades*a1/a2                            # Equação de conservação para encontrar a velocidade embaixo do carro
p2=p_atm+0.5*rho*velocidades**2 - 0.5*v2**2     # Equação de bernoulli para encontrar a pressão embaixo do carro
deltap=p_atm-p2
downforce=deltap*a_sup                          #calculo da Downforce decorrente da diferença de pressão
massa=downforce/9.8
cp=deltap/pdinamica
print(cp)

# Cálculo do drag e lift para cada velocidade
drags = drag(cd, pdinamica, af)
lifts = lift(ld, pdinamica,a_sup)
#drags_cooling= drag_cooling(rho,velocidades,a_cooling)
drags_frictions=drag_friction(velocidades,rho,mi,length,a_sup)

# Plotagem do gráfico
plt.figure(figsize=(10, 6))
plt.plot(velocidades, drags, label='Drag')
plt.plot(velocidades, lifts, label='Lift', linestyle='--')
#plt.plot(velocidades, drags_frictions, label='drag_frictions', linestyle=':')
#plt.plot(velocidades, downforce, label='downforce', linestyle='-.')

# Configurações do gráfico
plt.title('Gráfico de Drag e Lift em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
plt.show()

#Arrefecimento
vazao_ar= velocidades*a_cooling*rho #vazão massica das entradas de ar
deltaT=temp_f - temp_ar
calora=vazao_ar*c_ar*deltaT
