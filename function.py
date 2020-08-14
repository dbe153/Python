# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

def new_function(age,height):
    average=age**2+height
    return average

import numpy as np

x = np.arange(4) 
print("x =", x)
print("x + 5 =", x + 5) 
print("x - 5 =", x - 5) 
print("x * 2 =", x * 2)
print("x / 2 =", x / 2)
print("x // 2 =", x // 2) # floor division
print("-x = ", -x) 
print("x ** 2 = ", x ** 2)
print("x % 2 = ", x % 2)

from watson_developer_cloud import VisualRecognitionV3
visual_recognition = VisualRecognitionV3(
     version='2018-03-19',
     iam_apikey='p-d1c6242da82452e9a5fef253fa738f45fbe7b825'
     )
 