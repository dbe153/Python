# -*- coding: utf-8 -*-

'Spyder Editor'

import math
import numpy as np
import pandas as pd

#These values are inputs that can be changed
global Oil_rate
Oil_rate=3563 #bbl/d
global API
API=20
global Water_rate
Water_rate= 3563 #bbl/d
global Water_SG
Water_SG= 1.05
global Gas_rate
Gas_rate= 1000 #mmscf/d
global Gas_SG
Gas_SG= 0.58
global CO2
CO2=0
global H2S
H2S=0
global N2
N2=0
global TempF
TempF= 65
global TID
TID= 3.92 #Tubing ID inches

#These values are caluclated to correct for units
global GOR
GOR= Gas_rate/Oil_rate*1000
global Oil_SG
Oil_SG=141.5/(API+131.5)
global TempR
TempR=TempF+460
global Area
Area = (3.1416 * (TID / 12) ** 2) / 4

df = pd.DataFrame()
data = []
def Fluid_Calcs(P):
    global data
    global df
    global Gas_SG
    global Pb
    Pb = 15.7286 * ((GOR / Gas_SG) ** 0.7885) * ((10 ** (0.002 * TempF)) / (10 ** (0.0142 * API)))
    
    global Rs
    De_Ghetto_Heavy_Oil_A = Gas_SG * (1 + 0.1595 * API ** 0.4078 * TempF ** (-0.2466) * math.log10((P) / 114.7))
    De_Ghetto_Heavy_Oil_B = De_Ghetto_Heavy_Oil_A * (1 + 0.5912 * API * TempF * math.log10(P / 100) * 10 ** -4)
    Rs = ((De_Ghetto_Heavy_Oil_B * P ** 1.2057) / 56.434) * 10 ** (10.9267 * API / TempR)
    if Rs > GOR: 
        Rs = GOR
     
    global Rsw
    RswA = 8.15839 - 0.0612265 * TempF + 0.000191663 * TempF ** 2
    RswB = 0.0101021 - 0.0000744241 * TempF + 0.000000305553 * TempF ** 2 - 0.000000000294883 * TempF ** 3
    RswC = (-9.02505 + 0.130237 * TempF - 0.000853425 * TempF ** 2 + 0.00000234122 * TempF ** 3 - 0.00000000237049 * TempF ** 4) * 0.0000001
    Rsw = RswA + RswB * P + RswC * P ** 2
     
    global Bo
    De_Ghetto_A = GOR ** 0.755 * Gas_SG ** 0.25 * API ** -1.5 + 0.45 * TempF
    Bo = 0.98496 + 0.0001 * De_Ghetto_A ** 1.5
     
    global Bw
    Bw = 1 + 0.00012 * (TempF - 60) + 0.000001 * (TempF - 60) ** 2 - 0.00000333 * P
     
    global Z_Factor
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
    
    global Bg
    Bg=(5.04 * Z_Factor * (TempF + 459.62) / P) * 0.001
     
    global Live_Oil_Viscosity
    De_Ghetto_Dead_Oil = 10 ** (TempF ** (-1.163) * np.exp(6.9824 - 0.04658 * API)) - 1
    Live_Oil_Viscosity = (10.715 * (Bo + 100) ** (-0.515)) * De_Ghetto_Dead_Oil ** (5.44 * (Bo + 150) ** (-0.338))
    
    global Water_Viscosity
    Water_Viscosity=np.exp(1.003-0.01479*float(TempF)+0.00001982*float(TempF)**2)
     
    global Gas_Viscosity
    Gas_rho = 2.7 * Z_Factor * (P + 14.7) / (Z_Factor * (TempF + 460))
    Gas_Viscosity = (((9.4 + 0.02 * 28.966 * Gas_SG) * (TempF + 460) ** 1.5 / (209 + 19 * 28.966 * Gas_SG + (TempF + 460))) * 0.0001 * np.exp((3.5 + 986 / (TempF + 460) + 0.01 * 28.966 * Gas_SG) * (Gas_rho / 62.428) ** (2.4 - 0.2 * (3.5 + 986 / (TempF + 460) + 0.01 * 28.966 * Gas_SG))))
     
    global Oil_Density
    Gas_In_Solution = (((1.375 - 0.0002 * Rs) - (0.695 - 0.00008 * Rs)) / 34) * (API - 22) + (0.695 - 0.00008 * Rs)
    Oil_Density = (350 * Oil_SG + 0.0764 * Gas_In_Solution * Rs) / (5.615 * Bo)
     
    global Water_Density
    Water_Density = (62.4 * Water_SG) / Bw
    
    global Gas_Density
    Mol_Weight = Gas_SG * 28.97
    Gas_Density = P * Mol_Weight / Z_Factor / 10.73 / (TempR)
     
    global Oil_Interfacial_Tension
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
        Oil_Interfacial_Tension = 1 
    else:
        Oil_Interfacial_Tension = Oil_IT_A
    
    global Water_Interfacial_Tension
    Water_Interfacial_Tension = (1.42157E-21 * (P + 14.7) ** 6 - 4.32881E-17 * (P + 14.7) ** 5 + 4.75377E-13 * (P + 14.7) ** 4 - 0.00000000221001 * (P + 14.7) ** 3 + 0.00000423447 * (P + 14.7) ** 2 - 0.00870071 * (P + 14.7) + 52.0438) - (((1.42157E-21 * (P + 14.7) ** 6 - 4.32881E-17 * (P + 14.7) ** 5 + 4.75377E-13 * (P + 14.7) ** 4 - 0.00000000221001 * (P + 14.7) ** 3 + 0.00000423447 * (P + 14.7) ** 2 - 0.00870071 * (P + 14.7) + 52.0438) - (7.84314E-22 * (P + 14.7) ** 6 - 2.65743E-17 * (P + 14.7) ** 5 + 3.56966E-13 * (P + 14.7) ** 4 - 0.00000000243549 * (P + 14.7) ** 3 + 0.00000907989 * (P + 14.7) ** 2 - 0.0195933 * (P + 14.7) + 75.9343)) / (280 - 74)) * (280 - TempF)
    
    global VSL
    QOPT = Oil_rate * Bo * 5.614 / 86400
    QWPT = Water_rate * Bw * 5.614 / 86400
    QLPT = QOPT + QWPT
    VSL = QLPT / Area
     
    global VSG
    QGPT = (Oil_rate * (GOR - Rs) - Water_rate * Rsw) * Bg / 86400
    VSG = QGPT / Area
    if VSG < 0:
        VSG = 0
     
    global Oil_Res_Rate
    Oil_Res_Rate = Oil_rate * Bo
     
    global Water_Res_Rate
    Water_Res_Rate = Water_rate * Bw
    
    global Liquid_Res_Rate
    Liquid_Res_Rate = Oil_Res_Rate + Water_Res_Rate
     
    global Gas_Res_Rate
    if (Gas_rate - (Oil_rate * Rs + Water_rate * Rsw) * Bg) < 0:
        Gas_Res_Rate = 0
    else:
        Gas_Res_Rate = (Gas_rate - (Oil_rate * Rs + Water_rate * Rsw) * Bg)
     
    global Total_Res_Rate
    Total_Res_Rate = Liquid_Res_Rate + Gas_Res_Rate
     
    global Oil_in_Liquid
    '''Oil fraction in liquid at reservoir conditions'''
    Oil_in_Liquid = Oil_Res_Rate / Liquid_Res_Rate
     
    global Water_in_Liquid
    '''Water fraction in liquid at reservoir conditions'''
    Water_in_Liquid = Water_Res_Rate / Liquid_Res_Rate
     
    global Liquid_in_Total
    '''Liquid fraction in total fluid at reservoir conditions'''
    Liquid_in_Total = Liquid_Res_Rate / Total_Res_Rate
     
    global Gas_in_Total
    '''Gas fraction in total fluid at reservoir conditions'''
    Gas_in_Total =  Gas_Res_Rate / Total_Res_Rate
     
    global Liquid_Viscosity
    Liquid_Viscosity = Water_in_Liquid * Water_Viscosity + Oil_in_Liquid * Live_Oil_Viscosity
    
    global Liquid_Density
    Liquid_Density = Water_in_Liquid * Water_Density + Oil_in_Liquid * Oil_Density
    
    global Liquid_Interfacial_Tension
    Liquid_Interfacial_Tension = Water_in_Liquid * Water_Interfacial_Tension + Oil_in_Liquid * Oil_Interfacial_Tension
    
    global Fluid_Viscosity
    Fluid_Viscosity = Liquid_in_Total * Liquid_Viscosity + Gas_in_Total * Gas_Viscosity
    
    data.append([
        P,
        Pb,
        Rs,
        Rsw,
        Bo, 
        Bw,
        Z_Factor,
        Bg,
        Live_Oil_Viscosity,
        Water_Viscosity,
        Gas_Viscosity,
        Oil_Density,
        Water_Density,
        Gas_Density,
        Oil_Interfacial_Tension,
        Water_Interfacial_Tension,
        VSL,
        VSG,
        Oil_Res_Rate,
        Water_Res_Rate,
        Liquid_Res_Rate,
        Gas_Res_Rate,
        Total_Res_Rate,
        Oil_in_Liquid,
        Water_in_Liquid,
        Liquid_in_Total,
        Gas_in_Total,
        Liquid_Viscosity,
        Liquid_Density,
        Liquid_Interfacial_Tension,
        Fluid_Viscosity
        ])
    headers = ["P","Pb","Rs","Rsw","Bo","Bw","Z Factor","Bg","Live Oil Viscosity","Water Viscosity","Gas Viscosity","Oil Density","Water Density", "Gas Density", "Oil IT", "Water IT", "VSL","VSG","Qo_res","Qw_res","Ql_res","Qg_res","Qt_res","Qo_Ql","Qw_Ql","Ql_Qt","Qg_Qt","Liquid_cP","Liquid_Density","Liquid_IT","Fluid_cP"]
    df = pd.DataFrame(data)
    df.columns = headers
    return df
    print(df)


#write html to file
'''
html = df.to_html()
text_file = open("Fluid_Calcs.html", "w")
text_file.write(html)
text_file.close()
''' 



