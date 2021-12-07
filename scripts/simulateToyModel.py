# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 15:35:43 2021

@author: u0139894
"""

from pylab import *
import numpy as np
from gillespie_min import *
from BuildGillespieModel import Model
import matplotlib.pyplot as plt


dfs={}

for itr in range(20):
    M = Model('C:/Users/u0139894/Documents/GitHub/ChemSamp/files/pyruvate', 'test_')
    cpds = {i: Compound(M.mets[i][0], M.mets[i][1]) for i in M.mets}
    reacs = {i: Reaction(i, M.reacs[i]['rate'], {cpds[z]:M.reacs[i]['reactants'][z] for z in M.reacs[i]['reactants']}, {cpds[z]:M.reacs[i]['products'][z] for z in M.reacs[i]['products']}) for i in M.reacs}
    
    rsys = ReacSystem([reacs[i] for i in reacs])
    
    
    rsys.simul(200)
    
    df = pd.DataFrame.from_dict({rsys.cpdIds[i]:rsys.series.T[i] for i in range(len(rsys.cpdIds))})
    df.index = rsys.time
    #df.plot(title='k-leaps')
    
    #show()
    
    #df[["fru[out]", "tre[out]", "pyr", "pep", "f1p","fdp", "atp", "h", "dhap", "g3p", "nad", "13dpg"]].plot()
    dfs[itr] = df[['CO2']]

dd = dfs[0].copy()
for i in range(1,20):
    dd+=dfs[i]
dd=dd/20
dd.plot()