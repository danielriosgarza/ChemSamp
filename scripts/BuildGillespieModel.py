# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 08:48:06 2021

@author: u0139894
"""

import numpy as np
import os


class Model:
    def __init__(self, modelPath, initial):
        self.modelPath = modelPath
        self.initial = initial
        self.__getMets()
        self.__getReacs()
        
    
    def __getMets(self):
        mets = {}
        
        
        with open(os.path.join(self.modelPath, 'metabolites.txt')) as f:
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                if int(a[3]):
                    mets[a[0]] = (a[1], int(a[2]))
                else:
                    mets[self.initial+a[0]] = (a[1], int(a[2]))
        with open(os.path.join(self.modelPath, 'enzymes.txt')) as f:
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                mets[self.initial + a[0] + '_a'] = (a[1].upper() + '(free)', float(a[2]))
                mets[self.initial + a[0] + '_b'] = (a[1].upper() + '(conj)', 0)
        self.mets = mets
    
    def __makeString(self, st):
        sc = st.replace('"','')
        
        r = sc.split(',')
        rs = ''
        c = 0
        for i in r:
            a = i.split(':')
            
            
            if c==0:
                rs += str('(' + a[1] + ') ')
                
            else:
                rs+=' + ' + '(' + a[1] + ') '
            c+=1
            if a[0].replace(' ','') in self.mets:
                rs += self.mets[a[0].replace(' ','')][0]
            else:
                rs += self.mets[self.initial + a[0].replace(' ','')][0]
        return rs
    
    def __makeReacDict(self, st):
        d = {}
        
        sc = st.replace('"','')
        
        r = sc.split(',')
        
        
        for i,v in enumerate(r):
            a = v.split(':')
            met = a[0].replace(' ','')
            if self.initial + met in self.mets:
                met = self.initial + met
            
            d[met] = float(a[1])
        return d
        
        
    def __getReacs(self):
        reacs = {}
        with open(os.path.join(self.modelPath, 'reactions.txt')) as f:
            f.readline()
            for line in f:
                a = line.strip().split('\t')
                
                #conjugate
                reaction = a[0] + '_a'
                reactants = a[1] + ',' + a[3] + '_a:1'
                products = a[3] + '_b:1'
                
                reacs[reaction] = {'reactants':self.__makeReacDict(reactants), 'products': self.__makeReacDict(products), 'string':self.__makeString(reactants) + ' => ' + self.__makeString(products), 'rate':float(a[4])}
                
                #react
                reaction = a[0] + '_b'
                reactants = a[3] + '_b:1' 
                products = a[2] + ',' + a[3] + '_a:1'
                
                reacs[reaction] = {'reactants':self.__makeReacDict(reactants), 'products': self.__makeReacDict(products), 'string':self.__makeString(reactants) + ' => ' + self.__makeString(products), 'rate':float(a[4])}
                
                #rev conjugate
                reaction = a[0] + '_c'
                reactants = a[2] + ',' + a[3] + '_a:1'
                products = a[3] + '_b:1' 
                
                reacs[reaction] = {'reactants':self.__makeReacDict(reactants), 'products': self.__makeReacDict(products), 'string':self.__makeString(reactants) + ' => ' + self.__makeString(products), 'rate':float(a[5])}
                
                #rev react
                reaction = a[0] + '_d'
                reactants = a[3] + '_b:1'
                products =  a[1] + ',' + a[3] + '_a:1'
                
                reacs[reaction] = {'reactants':self.__makeReacDict(reactants), 'products': self.__makeReacDict(products), 'string':self.__makeString(reactants) + ' => ' + self.__makeString(products), 'rate':float(a[5])}

                    
                    
                
                
        self.reacs = reacs
        

#M = Model('C:/Users/u0139894/Documents/BH_EnergyMetabolism/files/gillespie/autotrophic', 'aut_')

# with open('C:/Users/u0139894/Documents/BH_EnergyMetabolism/files/gillespie/reactions.txt') as f:
#     with open('C:/Users/u0139894/Documents/BH_EnergyMetabolism/files/gillespie/reactions2.txt', 'w') as f2:
#         f2.write(f.readline())
#         for line in f:
#             a=line.strip().split('\t')
#             f2.write(line.strip() + '\t' + M.reacs[a[0]] + '\n')
            
        