# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 20:31:03 2020

@author: Lennard
"""


import numpy as np
from loadcurves import loadcurves as load

#material Parameters
E=69*10**9
mue=0.33
H=2000*10**6
M=8
a=1.24
b=1.02
h=1.15
G=25.5*10**9
sigyield11=1

mat = np.array([E,mue,H,M,a,b,h,G,sigyield11])

#choosing of a load curve
#1 linear increasing stress(zero -> max)
#2 cyclic stress state (0 -> max -> minus max -> max -> minus max -> max)
    
flag1 = 2

#amplitude of eps11 and eps22
eps11_max = 1
eps22_max = 3

if flag1 == 1:
    loadcurve=load(flag1,eps11_max,eps22_max)
elif flag1 == 2:
    loadcurve=load(flag1,eps11_max,eps22_max)
else:
    print('Choose a valid value for flag1')

    

























    