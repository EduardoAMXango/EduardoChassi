import numpy as np
import matplotlib.pyplot as plt
import math

#Parâmetros físicos do ar:
p_atm= 101325                       #pressâo (N/m^2)
rho= 1.184                          #densidade(kg/m^3)
mi= 1.849*10**-5                    #Viscosidade dinâmica (kg/m*s)
ni= (mi/rho)                        #Viscosidade Cinematica (m²/s)
c_ar= 1007                          #Calor especifico (J/(kg*K))
#parâmetros do carro
length= 2                           #Comprimento do carro (m^2)
af= 1.5                             #Área Frontal do carro (m^2)
a_sup= 2                            #Área de Superfície do carro (m^2)
cd_f = 0.05                         #Coeficiente de arrasto por atrito do carro
cd_p = 0.45                         #Coeficiente de arrasto por pressão do carro
ld= -0.3                            #Coeficiente de lift do carro
a1= 0.25                            #Área de entrada do ar embaixo do carro (m^2)
a2= 0.20                            #Área embaixo do carro (m^2)
#parâmetros da asa
chord=0.25                          #Comprimento da asa (m)
span=1                              #Largura da asa (m)
thickness=0.05                      #Expessura máxima da asa (m)
alpha=math.radians(-3.75)               #Ângulo de incidencia do vento com a asa (radianos)
a_wing=span*chord                   #Área superior da asa (m^2)
aspect_ratio=(span**2)/a_wing       #Comparação entre a largura e o comprimento da asa
#parâmetros para a troca termica
deltaT=20                            #Variação de temperatura (K)
a_ent=0.05                           #Área de admissão do vento para arrefecimento (m^2)
a_res = 1                            #Área de resfriamento (m^2)
h= 300                               #Coeficiente de transferência de calor (W/m2*K)

#Equações Aerodinâmicas: 
def aerodynamic_forces(coeficient, pressao_dinamica,area):       # Função para o cálculo de arrasto e lift
    return coeficient * pressao_dinamica * area

def numero_Reynolds(comprimento,velocidade,ni):         # Função para encontrar o numero de reynolds
    return (comprimento*velocidade)/ni

def bernoulli(v1,p1,a1,a2,rho):
    v2=v1*a1/a2                            # Equação de conservação para encontrar a velocidade embaixo do carro
    p2=p1+0.5*rho*v1**2 - 0.5*v2**2        # Equação de bernoulli para encontrar a pressão embaixo do carro
    deltap=p1-p2                           # Diferença de pressão entre a parte de baixo e a parte de cima do carro
    return deltap

def lift_wing(alpha,aspect_ratio):
    cl_alpha=(2*math.pi)/(1+(2/aspect_ratio))
    cl_wing=cl_alpha*alpha
    return cl_wing

#Soma dos coeficientes de arrasto
cd = cd_f+cd_p

# Velocidades de 0 a 60 m/s
velocidades = np.linspace(0.1,60,50)
# Cálculo dos coeficientes das asas
cl_w=lift_wing(alpha,aspect_ratio)                  #Coeficiente de lift da asa
cd_induced=(1/(math.pi*aspect_ratio))*cl_w**2       #Coeficiente de arrasto induzido da asa
cd_w=cd_induced+0.05
print(f"coeficiente de lift da asa: {cl_w:0.2f} , coeficiente de arrasto induzido: {cd_induced:0.2f}")
for v in velocidades:
    # Cálculo da pressão dinâmica para cada velocidade
    pdinamica = 0.5 * rho  * (velocidades ** 2)

    # Cálculo do drag e lift para cada velocidade
    drags = aerodynamic_forces(cd_p, pdinamica, af) + aerodynamic_forces(cd_f, pdinamica, a_sup)                                    #Arrasto 
    downforces = aerodynamic_forces(ld, pdinamica,a_sup)                            #Downforce 

    lifts_wing = aerodynamic_forces(cl_w, pdinamica, a_wing)                        #Lift da asa
    drags_wing = aerodynamic_forces(cd_w, pdinamica, a_wing)                  #Arrasto induzido pela asa

    #Calculo de Downforce por bernoulli
    downforce_bernoulli= (-1)*bernoulli(velocidades,p_atm,a1,a2,rho)*a_sup    #Conversão da pressão para força por meio da relação pressão=força/área --- força=pressão*área

# Plotagem do gráfico 1
plt.figure(figsize=(10, 6))
plt.subplot(1,3,1)
plt.plot(velocidades, drags, label='Drag')
# Configurações do gráfico
plt.title('Gráfico de Arrasto em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)

# Plotagem do gráfico 2
plt.subplot(1,3,2)
plt.plot(velocidades, downforces, label='Downforce', linestyle='--')
plt.plot(velocidades, downforce_bernoulli, label='Downforce por bernoulli', linestyle='-.')
# Configurações do gráfico2
plt.title('Gráfico de Downforce em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
# Plotagem do gráfico 3
plt.subplot(1,3,3)
plt.plot(velocidades, lifts_wing, label='Lift', linestyle='--')
plt.plot(velocidades, drags_wing, label='Arrasto', linestyle='-.')
# Configurações do gráfico 3
plt.title('Gráfico de forças da asa em função da Velocidade')
plt.xlabel('Velocidade (m/s)')
plt.ylabel('Força (N)')
plt.legend()
plt.grid(True)
plt.show()

