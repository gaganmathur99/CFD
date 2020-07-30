# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:24:48 2019

@author: Arpit Mathur
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import time

start=time.time();
T1=100
T2=1000;
Lx=1;
Ly=1;
div=1000;
x=np.linspace(0,Lx,div);
y=np.linspace(0,Ly,div);

theta=np.zeros([div,div]);
T=np.zeros([div,div]);

for j in range (0,div):
    i=0;
    T[i][j]=T1;
for j in range (0,div-1):
    i=div-1;
    T[i][j]=T1;
for i in range (0,div):
    j=0;
    T[i][j]=T1;
for i in range (0,div-1):
    j=div-1;
    T[i][j]=T2;


for n in range (1,20):
    for i in range (1,div):
        for j in range (1,div):
            theta[i][j]=(2/math.pi)*((((-1)**n)+1)/n)*(math.sin(n*math.pi*x[i]/Lx))*math.asinh(n*math.pi*y[j]/Lx)/math.asinh(n*math.pi*Ly/Lx);
            T[i][j]=theta[i][j]*(T2-T1)+T1;
    print(n)
plt.contourf(T)