# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 03:45:20 2021

@author: danie
"""

from pylab import *
import numpy as np

def gPh(pH, n1, n2):
    pHmin, pHmax = 6.0, 10
    return (pH**n1)/(((pHmin**n1)-1) + (pH**n1)) * (pHmax**n2)/(((pHmax**n2)-1) + (pH**n2))

n1=25
n2 = 0

gPlot = np.array([gPh(i, n1, n2) for i in np.linspace(4,11,1000)])
