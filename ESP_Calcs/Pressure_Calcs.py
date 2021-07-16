#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 16:21:57 2021

@author: dellexson
"""
#import os
import math
#import numpy as np
import pandas as pd
import Fluid_Calcs_Formulas
#import Deviation_Surveys

#print(Deviation_Surveys)

Fluid_Calcs_Formulas.Fluid_Calcs(200)



if VSL < 0.1:
    VSL = 0.1
VM = VSL + VSG
if Radians < 0.01:
    Radians = 0.01

WOR = Q_Water / Q_Oil
Produced_GLR = GLR - Rso - Rsw
if Produced_GLR < 0.01:
    Produced_GLR = 0
Fluid_Mass = (Oil_SG * 350 * (1 / (1 + WOR))) + (Water_SG * 350 * (WOR / (1 + WOR))) + (0.0764 * Produced_GLR * Gas_SG)
Mass_Flow_Rate = Fluid_Mass * (Q_Oil + Q_Water)

NL = 0.15726 * Liquid_Visc * (1 / (Liquid_Density * stl ** 3)) ** 0.25
CNL = 10 ** (-2.69851 + 0.15840954 * (WorksheetFunction.Log10(NL) + 3) + -0.55099756 * (math.log10(NL) + 3) ^ 2 + 0.54784917 * (math.log10(NL) + 3) ^ 3 + -0.12194578 * (WorksheetFunction.Log10(NL) + 3) ^ 4)
NLV = (1.938 * VSL * ((Liquid_Density / stl)) ^ (0.25))
if (1.938 * VSG * ((Liquid_Density / stl)) ** (0.25)) < 0.01:
    NGV = 0.01 
else:
    NGV = (1.938 * VSG * ((Liquid_Density / stl)) ** (0.25))
Nd = (120.872 * (TID / 12) * ((Liquid_Density / stl)) ^ (0.5))
CNL1 = NLV / NGV ^ 0.575 * (segment_psi / 14.7) ^ 0.1 * CNL / Nd

def Hagedorn_Brown_Flow_Pattern():
    global TID
    global Area
    global VSL
    global VSG
    if 1.071 - (0.2218 * (VSL + VSG) ** 2) / (TID / 12) < 0.13:
        A_HB = 0.13 
    else:
        A_HB = 1.071 - (0.2218 * (VSL + VSG) ** 2) / (TID / 12)
        B_HB = VSG / (VSL + VSG)
        if B_HB - A_HB < 0:
            Hagedorn_Brown_Flow_Pattern = "Bubble" 
        else:
            Hagedorn_Brown_Flow_Pattern = "Slug"
        if B_HB - A_HB <= 0:
            HG = (1 / 2) * (1 + (Liquid_Q + Q_Gas_BPD) * 0.0000649836034 / (0.8 * Area) - (((1 + (Liquid_Q + Q_Gas_BPD) * 0.0000649836034 / (0.8 * Area)) ** 2 - 4 * Q_Gas_BPD * 0.0000649836034 / (0.8 * Area)) ** 0.5))
            HL_G = 1 - HG
            Bre = 1488 * Liquid_Density * (TID / 12) * (Liquid_Q * 0.0000649836034 / Area) / Liquid_Visc
            if Bre >= 2000:
                BF = (-2 * math.log10((ed / (TID / 12)) / 3.7 - 5.02 / Bre * math.log10((ed / (TID / 12)) / 3.7 + 13 / Bre))) ** (-2) 
            else:
                BF = 64 / Bre
            Elev = 1 / 144 * (32.1522 * Liquid_Density * Sin(Radians) / 32.174)
            Fric = 1 / 144 * ((BF * Liquid_Density * VSL ^ 2) / (2 * 32.174 * (TID / 12)))
            
            dpdz = Elev + Fric
        
        else:
            Holdup_C = (NLV / NGV ** 0.575) * ((segment_psi / 14.7) ^ 0.1) * (CNL / Nd)
            HL_Holdup_C = -0.10306578 + 0.617774 * (math.log10(Holdup_C) + 6) + -0.632946 * (math.log10(Holdup_C) + 6) ** 2 + 0.29598 * (math.log10(Holdup_C) + 6) ** 3 + -0.0401 * (math.log10(Holdup_C) + 6) ** 4
            Abscissa = NGV * NL ** 0.38 / Nd ** 2.14
            Y = 0.91162574 - 4.82175636 * Abscissa + 1232.25036621 * Abscissa ^ 2 - 22253.57617 * Abscissa ^ 3 + 116174.28125 * Abscissa ^ 4
            HL_HB = HL_Holdup_C / Y
            Nre = (2.2 * 10 ^ -2) * Mass_Flow_Rate / ((TID / 12) * Liquid_Visc ** HL_HB * gas_Visc ** (1 - HL_HB))
            ed_HB = (ed * 12) / TID
            Fanning_F = 1 / (-4 * math.log10(ed_HB / 3.7065 - 5.0452 / Nre * math.log10(ed_HB ^ 1.1098 / 2.8257 + (7.149 / Nre) ^ 0.8981))) ^ 2
            p_g = 28.97 * Gas_SG * segment_psi / Z_HB / 10.73 / (Temp_R)
            Ave_P = (HL_HB * Liquid_Density + (1 - HL_HB) * p_g) * Sin(Radians)
            
            dpdz = 1 / 144 * (Ave_P + Fanning_F * Mass_Flow_Rate ** 2 / 7.413 / 10000000000 / (TID / 12) ** 5 / Ave_P)



for i in range(200,1600,100):
    table = Fluid_Calcs_Formulas.Fluid_Calcs(i)
    i + 200
print(table)
