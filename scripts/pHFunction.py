# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 03:45:20 2021

@author: danie
"""

from pylab import *
import numpy as np

def gPh(pH, n1, n2, pHmin, pHmax):
    return (pH**n1)/(((pHmin**n1)-1) + (pH**n1)) * (pHmax**n2)/(((pHmax**n2)-1) + (pH**n2))

n1=12
n2 = 5
pHmin = 5
pHmax = 9

gPh = np.array([gPh(i, n1, n2, pHmin, pHmax) for i in np.linspace(4,10,1000)])
