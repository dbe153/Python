# -*- coding: utf-8 -*-

'Spyder Editor'

import math
import numpy as np
import pandas as pd
from IPython.display import display

#These values are inputs that can be changed
Oil_rate=3563 #bbl/d
API=20
Water_rate= 3563 #bbl/d
Water_SG= 1.05
Gas_rate= 1000 #mmscf/d
Gas_SG= 0.58
CO2=0
H2S=0
N2=0
TempF= 65
TID= 3.92 #Tubing ID inches

#These values are caluclated to correct for units
GOR= Gas_rate/Oil_rate*1000
Oil_SG=141.5/(API+131.5)
TempR=TempF+460
Area = (3.1416 * (TID / 12) ** 2) / 4

data = []
for P in (200,372,544,1500):

    def De_Ghetto_Pb(TempF,API,Gas_SG,GOR):
         '''Bubble Point Pressure(psi)'''
         Oil_Pb=15.7286 * ((GOR / Gas_SG) ** 0.7885) * ((10 ** (0.002 * TempF)) / (10 ** (0.0142 * API)))
         return Oil_Pb
    
    def De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR):
         De_Ghetto_Heavy_Oil_A = Gas_SG * (1 + 0.1595 * API ** 0.4078 * TempF ** (-0.2466) * math.log10((P) / 114.7))
         De_Ghetto_Heavy_Oil_B = De_Ghetto_Heavy_Oil_A * (1 + 0.5912 * API * TempF * math.log10(P / 100) * 10 ** -4)
         De_Ghetto_Heavy_Oil_Rs = ((De_Ghetto_Heavy_Oil_B * P ** 1.2057) / 56.434) * 10 ** (10.9267 * API / TempR)
     
         if De_Ghetto_Heavy_Oil_Rs > GOR: 
             De_Ghetto_Heavy_Oil_Rs = GOR
         return De_Ghetto_Heavy_Oil_Rs
     
    def Water_Rsw(TempF,P):
         RswA = 8.15839 - 0.0612265 * TempF + 0.000191663 * TempF ** 2
         RswB = 0.0101021 - 0.0000744241 * TempF + 0.000000305553 * TempF ** 2 - 0.000000000294883 * TempF ** 3
         RswC = (-9.02505 + 0.130237 * TempF - 0.000853425 * TempF ** 2 + 0.00000234122 * TempF ** 3 - 0.00000000237049 * TempF ** 4) * 0.0000001
         Rsw = RswA + RswB * P + RswC * P ** 2
         return Rsw
     
    def De_Ghetto_Bo(TempF,Gas_SG,GOR,API):
         De_Ghetto_A = GOR ** 0.755 * Gas_SG ** 0.25 * API ** -1.5 + 0.45 * TempF
         Bo = 0.98496 + 0.0001 * De_Ghetto_A ** 1.5
         return Bo
     
    def Gould_Bw(TempF,P):
         Bw = 1 + 0.00012 * (TempF - 60) + 0.000001 * (TempF - 60) ** 2 - 0.00000333 * P
         return Bw
     
    def Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P):
         if Gas_SG > 1.25:
             Gas_SG = 1.25
         Pseudocritical_pressure = 678 - 50 * (Gas_SG - 0.5) - 206.7 * N2 + 440 * CO2 + 606.7 * H2S
         Pseudocritical_temperature = 326 + 315.7 * (Gas_SG - 0.5) - 240 * N2 - 83.3 * CO2 + 133.3 * H2S
         Pseudo_reduced_pressure = P / Pseudocritical_pressure
         pseudo_reduced_temperature = (TempR) / Pseudocritical_temperature
         a = 1.39 * (pseudo_reduced_temperature - 0.92) ** 0.5 - 0.36 * pseudo_reduced_temperature - 0.101
         b = (0.62 - 0.23 * pseudo_reduced_temperature) * Pseudo_reduced_pressure + (0.066 / (pseudo_reduced_temperature - 0.86) - 0.037) * Pseudo_reduced_pressure ** 2 + 0.32 / 10 ** (9 * (pseudo_reduced_temperature - 1)) * Pseudo_reduced_pressure ** 6
         c = 0.132 - 0.32 * math.log10(pseudo_reduced_temperature)
         d = 10 ** (0.3106 - 0.49 * pseudo_reduced_temperature + 0.1824 * pseudo_reduced_temperature ** 2)
         Z_Factor = a + (1 - a) / np.exp(b) + c * Pseudo_reduced_pressure ** d
         return Z_Factor
    
    def Standing_Bg(TempF,Beggs_Brill_Z_Factor,P):
         Bg=(5.04 * Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P) * (TempF + 459.62) / P) * 0.001
         return Bg
     
    def De_Ghetto_Viscosity(TempF,API,Gas_SG,GOR,P,De_Ghetto_Bo):
         De_Ghetto_Dead_Oil = 10 ** (TempF ** (-1.163) * np.exp(6.9824 - 0.04658 * API)) - 1
         De_Ghetto_cP = (10.715 * (De_Ghetto_Bo(TempF,Gas_SG,GOR,API) + 100) ** (-0.515)) * De_Ghetto_Dead_Oil ** (5.44 * (De_Ghetto_Bo(TempF,Gas_SG,GOR,API) + 150) ** (-0.338))
         return De_Ghetto_cP
    
    def Van_Wingen_Water_Viscosity(TempF):
         Water_Viscosity=np.exp(1.003-0.01479*float(TempF)+0.00001982*float(TempF)**2)
         return Water_Viscosity
     
    def Lee_Gas_Viscosity(TempF,P,Beggs_Brill_Z_Factor):
         Gas_rho = 2.7 * Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P) * (P + 14.7) / (Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P) * (TempF + 460))
         Gas_Viscosity = (((9.4 + 0.02 * 28.966 * Gas_SG) * (TempF + 460) ** 1.5 / (209 + 19 * 28.966 * Gas_SG + (TempF + 460))) * 0.0001 * np.exp((3.5 + 986 / (TempF + 460) + 0.01 * 28.966 * Gas_SG) * (Gas_rho / 62.428) ** (2.4 - 0.2 * (3.5 + 986 / (TempF + 460) + 0.01 * 28.966 * Gas_SG))))
         return Gas_Viscosity
     
    def Oil_Density(De_Ghetto_Rs,API,De_Ghetto_Bo,Oil_SG):
         Gas_In_Solution = (((1.375 - 0.0002 * De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR)) - (0.695 - 0.00008 * De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR))) / 34) * (API - 22) + (0.695 - 0.00008 * De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR))
         oil_rho = (350 * Oil_SG + 0.0764 * Gas_In_Solution * De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR)) / (5.615 * De_Ghetto_Bo(TempF,Gas_SG,GOR,API))
         return oil_rho
     
    def Water_Density(Water_SG,Gould_Bw):
         Water_density = (62.4 * Water_SG) / Gould_Bw(TempF,P)
         return Water_density
     
    def Gas_Density(Gas_SG,TempR,Beggs_Brill_Z_Factor):
         Mol_Weight = Gas_SG * 28.97
         Gas_density = P * Mol_Weight / Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P) / 10.73 / (TempR)
         return Gas_density
     
    def Oil_Interfacial_Tension(API, TempF, P):
         ST68 = 39 - 0.2571 * API
         ST100 = 37.5 - 0.2571 * API
         tst = TempF
         if TempF > 100:
             tst = 100
         if TempF < 68:
             Oil_IT_A = ST68 
         else:
             Oil_IT_A = (68 - (((tst - 68) * (ST68 - ST100)) / 32)) * (1 - (0.024 * (P + 14.5) ** 0.45))
         if Oil_IT_A < 1:
             Oil_IT = 1 
         else:
             Oil_IT = Oil_IT_A
             
         return Oil_IT
     
    def Water_Interfacial_Tension(TempF,P,):
         Water_IT = (1.42157E-21 * (P + 14.7) ** 6 - 4.32881E-17 * (P + 14.7) ** 5 + 4.75377E-13 * (P + 14.7) ** 4 - 0.00000000221001 * (P + 14.7) ** 3 + 0.00000423447 * (P + 14.7) ** 2 - 0.00870071 * (P + 14.7) + 52.0438) - (((1.42157E-21 * (P + 14.7) ** 6 - 4.32881E-17 * (P + 14.7) ** 5 + 4.75377E-13 * (P + 14.7) ** 4 - 0.00000000221001 * (P + 14.7) ** 3 + 0.00000423447 * (P + 14.7) ** 2 - 0.00870071 * (P + 14.7) + 52.0438) - (7.84314E-22 * (P + 14.7) ** 6 - 2.65743E-17 * (P + 14.7) ** 5 + 3.56966E-13 * (P + 14.7) ** 4 - 0.00000000243549 * (P + 14.7) ** 3 + 0.00000907989 * (P + 14.7) ** 2 - 0.0195933 * (P + 14.7) + 75.9343)) / (280 - 74)) * (280 - TempF)
         return Water_IT
     
    def Superficial_Liquid(TID,Oil_rate,De_Ghetto_Bo,Water_rate,Gould_Bw,Area):
         QOPT = Oil_rate * De_Ghetto_Bo(TempF,Gas_SG,GOR,API) * 5.614 / 86400
         QWPT = Water_rate * Gould_Bw(TempF,P) * 5.614 / 86400
         QLPT = QOPT + QWPT
         VSL = QLPT / Area
         return VSL
     
    def Superficial_Gas(TID,Oil_rate,De_Ghetto_Rs,Water_rate,Water_Rsw,Standing_Bg,Area):
         QGPT = (Oil_rate * (GOR - De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR)) - Water_rate * Water_Rsw(TempF,P)) * Standing_Bg(TempF,Beggs_Brill_Z_Factor,P) / 86400
         VSG = QGPT / Area
         if VSG < 0:
             VSG = 0
         return VSG
     
    def Oil_Res_Rate(Oil_rate,De_Ghetto_Bo):
          Qo_res = Oil_rate * De_Ghetto_Bo(TempF,Gas_SG,GOR,API)
          return Qo_res
    #print("OIl_res",Oil_res_rate(Oil_rate,De_Ghetto_Bo(TempF,Gas_SG,GOR,API)))
     
    def Water_Res_Rate(Water_rate,Gould_Bw):
         Qw_res = Water_rate * Gould_Bw(TempF,P)
         return Qw_res
     
    def Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate):
         Ql_res = Oil_Res_Rate(Oil_rate,De_Ghetto_Bo) + Water_Res_Rate(Water_rate,Gould_Bw)
         return Ql_res
     
    def Gas_Res_Rate(Gas_rate,Oil_rate,Water_rate,De_Ghetto_Rs,Water_Rsw,Standing_Bg):
         if (Gas_rate - (Oil_rate * De_Ghetto_Rs(TempF, API, Gas_SG, GOR, P, TempR) + Water_rate * Water_Rsw(TempF,P)) * Standing_Bg(TempF,Beggs_Brill_Z_Factor,P)) < 0:
             Qg_res = 0
         else:
              (Gas_rate - (Oil_rate * De_Ghetto_Rs(TempF, API, Gas_SG, GOR, P, TempR) + Water_rate * Water_Rsw(TempF,P)) * Standing_Bg(TempF,Beggs_Brill_Z_Factor,P))
         return Qg_res
     
    def Total_Res_Rate(Liquid_Res_Rate,Gas_Res_Rate):
         Qt_res = Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate) + Gas_Res_Rate(Gas_rate,Oil_rate,Water_rate,De_Ghetto_Rs,Water_Rsw,Standing_Bg)
         return Qt_res
     
    def Oil_in_Liquid(Oil_Res_Rate,Liquid_Res_Rate):
         '''Oil fraction in liquid at reservoir conditions'''
         Qo_Ql = Oil_Res_Rate(Oil_rate,De_Ghetto_Bo) / Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate)
         return Qo_Ql
     
    def Water_in_Liquid(Water_Res_Rate,Liquid_Res_Rate):
         '''Water fraction in liquid at reservoir conditions'''
         Qw_Ql = Water_Res_Rate(Water_rate,Gould_Bw) / Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate)
         return Qw_Ql
     
    def Liquid_in_Total(Liquid_Res_Rate,Total_Res_Rate):
         '''Liquid fraction in total fluid at reservoir conditions'''
         Ql_Qt = Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate) / Total_Res_Rate(Liquid_Res_Rate,Gas_Res_Rate)
         return Ql_Qt
     
    def Gas_in_Total(Gas_Res_Rate,Total_Res_Rate):
         '''Gas fraction in total fluid at reservoir conditions'''
         Qg_Qt =  Gas_Res_Rate(Gas_rate,Oil_rate,Water_rate,De_Ghetto_Rs,Water_Rsw,Standing_Bg) / Total_Res_Rate(Liquid_Res_Rate,Gas_Res_Rate)
         return Qg_Qt
     
    def Liquid_Viscosity(Water_in_Liquid,Oil_in_Liquid,De_Ghetto_Viscosity,Van_Wingen_Water_Viscosity):
        Liquid_cP = Water_in_Liquid(Water_Res_Rate,Liquid_Res_Rate) * Van_Wingen_Water_Viscosity(TempF) + Oil_in_Liquid(Oil_Res_Rate,Liquid_Res_Rate) * De_Ghetto_Viscosity(TempF,API,Gas_SG,GOR,P,De_Ghetto_Bo)
        return Liquid_cP
    
    def Liquid_Density(Water_in_Liquid,Oil_in_Liquid,Oil_Density,Water_Density):
        Liquid_density = Water_in_Liquid(Water_Res_Rate,Liquid_Res_Rate) * Water_Density(Water_SG,Gould_Bw) + Oil_in_Liquid(Oil_Res_Rate,Liquid_Res_Rate) * Oil_Density(De_Ghetto_Rs,API,De_Ghetto_Bo,Oil_SG)
        return Liquid_density
    
    def Liquid_Interfacial_Tension(Water_in_Liquid,Oil_in_Liquid,Oil_Interfacial_Tension,Water_Interfacial_Tension):
        Liquid_IT = Water_in_Liquid(Water_Res_Rate,Liquid_Res_Rate) * Water_Interfacial_Tension(TempF,P,) + Oil_in_Liquid(Oil_Res_Rate,Liquid_Res_Rate) * Oil_Interfacial_Tension(API, TempF, P)
        return Liquid_IT
    
    def Fluid_Viscosity(Liquid_in_Total,Gas_in_Total,Liquid_Viscosity,Lee_Gas_Viscosity):
        Fluid_cP = Liquid_in_Total(Liquid_Res_Rate,Total_Res_Rate) * Liquid_Viscosity(Water_in_Liquid,Oil_in_Liquid,De_Ghetto_Viscosity,Van_Wingen_Water_Viscosity) + Gas_in_Total(Gas_Res_Rate,Total_Res_Rate) * Lee_Gas_Viscosity(TempF,P,Beggs_Brill_Z_Factor)
        return Fluid_cP
    
    def Fluid_Density():
        return Fluid_density
    
    def VM:
        return Mixture_Velocity
    
    def Ens:
        return 
    
    def Liquid_Velocity_Number():
        LVN = Superficial_Liquid(TID,Oil_rate,De_Ghetto_Bo,Water_rate,Gould_Bw,Area) * 0.3048 * (Liquid_Density(Water_in_Liquid,Oil_in_Liquid,Oil_Density,Water_Density) * 16.018463306 / Liquid_Interfacial_Tension(Water_in_Liquid,Oil_in_Liquid,Oil_Interfacial_Tension,Water_Interfacial_Tension) / 0.001 / 9.81)**0.25
        return LVN
    
    data.append([P,De_Ghetto_Pb(TempF,API,Gas_SG,GOR),De_Ghetto_Rs(TempF,API,Gas_SG,GOR,P,TempR), Water_Rsw(TempF,P),De_Ghetto_Bo(TempF,Gas_SG,GOR,API), Gould_Bw(TempF,P),\
         Beggs_Brill_Z_Factor(TempF,TempR,Gas_SG,CO2,N2,H2S,P),Standing_Bg(TempF,Beggs_Brill_Z_Factor,P),De_Ghetto_Viscosity(TempF,API,Gas_SG,GOR,P,De_Ghetto_Bo),\
             Van_Wingen_Water_Viscosity(TempF),Lee_Gas_Viscosity(TempF,P,Beggs_Brill_Z_Factor),Oil_Density(De_Ghetto_Rs,API,De_Ghetto_Bo,Oil_SG),Water_Density(Water_SG,Gould_Bw),\
                 Gas_Density(Gas_SG,TempR,Beggs_Brill_Z_Factor),Oil_Interfacial_Tension(API,TempF,P),Water_Interfacial_Tension(TempF,P),Superficial_Liquid(TID,Oil_rate,De_Ghetto_Bo,\
                     Water_rate,Gould_Bw,Area),Superficial_Gas(TID,Oil_rate,De_Ghetto_Rs,Water_rate,Water_Rsw,Standing_Bg,Area),Oil_Res_Rate(Oil_rate,De_Ghetto_Bo),\
                         Water_Res_Rate(Water_rate,Gould_Bw),Liquid_Res_Rate(Oil_Res_Rate,Water_Res_Rate),Gas_Res_Rate(Gas_rate,Oil_rate,Water_rate,De_Ghetto_Rs,Water_Rsw,Standing_Bg),\
                              Total_Res_Rate(Liquid_Res_Rate,Gas_Res_Rate),Oil_in_Liquid(Oil_Res_Rate,Liquid_Res_Rate),Water_in_Liquid(Water_Res_Rate,Liquid_Res_Rate),\
                                   Liquid_in_Total(Liquid_Res_Rate,Total_Res_Rate),Gas_in_Total(Gas_Res_Rate,Total_Res_Rate),Liquid_Viscosity(Water_in_Liquid,Oil_in_Liquid,De_Ghetto_Viscosity,Van_Wingen_Water_Viscosity),\
                                   Liquid_Density(Water_in_Liquid,Oil_in_Liquid,Oil_Density,Water_Density),Liquid_Interfacial_Tension(Water_in_Liquid,Oil_in_Liquid,Oil_Interfacial_Tension,Water_Interfacial_Tension),\
                                   Fluid_Viscosity(Liquid_in_Total,Gas_in_Total,Liquid_Viscosity,Lee_Gas_Viscosity)])

headers = ["P","Pb","Rs","Rsw","Bo","Bw","Z Factor","Bg","Live Oil Viscosity","Water Viscosity","Gas Viscosity","Oil Density","Water Density", "Gas Density", "Oil IT", "Water IT", "VSL","VSG","Qo_res","Qw_res","Ql_res","Qg_res","Qt_res","Qo_Ql","Qw_Ql","Ql_Qt","Qg_Qt","Liquid_cP","Liquid_Density","Liquid_IT","Fluid_cP"]
    
df = pd.DataFrame(data)
df.columns = headers
display(df)


#write html to file
'''
html = df.to_html()
text_file = open("Fluid_Calcs.html", "w")
text_file.write(html)
text_file.close()
''' 



