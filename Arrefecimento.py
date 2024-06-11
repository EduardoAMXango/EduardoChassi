import numpy as np
import matplotlib.pyplot as plt
import math

class Arrefecimento():
    def __init__(self,area_resfriamento,area_entrada,densidade_ar,coeficiente_transferencia,calor_especifico_ar):
        self.a_res= area_resfriamento
        self.a_ent= area_entrada
        self.rho= densidade_ar
        self.h= coeficiente_transferencia
        self.c= calor_especifico_ar
    def dissipação(self,deltaT):
        Qconvecção= self.h*self.a_res*deltaT
        return Qconvecção
    def calorimetria(self,velocidade,deltaT):
        calor=self.c*deltaT*(self.rho*self.a_ent*velocidade)
        return calor

#Parâmetros físicos do ar:
rho= 1.184                       #densidade(kg/m^3)
c_ar= 1007                       #Calor especifico (J/(kg*K))
#parâmetros para a troca termica
deltaT=20                        #Variação de temperatura (K)
a_ent=0.05                       #Área de admissão do vento para arrefecimento (m^2)
a_res = 1                        #Área de resfriamento (m^2)
h= 300                           #Coeficiente de transferência de calor (W/m2*K)
v=10
bateria=Arrefecimento(a_res,a_ent,rho,h,c_ar)
print(f"{bateria.dissipação(deltaT):0.2f}W")
print(f"{bateria.calorimetria(v,deltaT):0.2f}W")