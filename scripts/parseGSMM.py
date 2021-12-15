# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 10:38:16 2021

@author: u0139894
"""

import cobra
import numpy as np

######template model#########


def getReacProd(reactionObj, el='reactants'):
    string = ''
    d_el = reactionObj.reactants
    
    
    if el=='products':
        d_el = reactionObj.products
    
    if len(d_el)==0:
        if el=='products':
            d_el = reactionObj.reactants
            
        elif el=='reactants':
            d_el = reactionObj.products
            
            
        
    for i in d_el:
        string += i.id
        string += ':'
        string += str(abs(reactionObj.metabolites[i]))
        string += ','
    return string[0:-1]


def extractReactions(templateModel, reactionList, modelRoot):
    
    metabolites = {}
    reactions = {}
    enzymes = {}
    
    for reaction in reactionList:
        
        reac = model.reactions.get_by_id(reaction)
        for metabolite in reac.metabolites:
            metabolites[metabolite.id] = {'id':metabolite.id, 'abbrev':metabolite.notes['abbrev'], 'n':0, 'exchange': int('e' in metabolite.compartment)}
            
            
        enzymes[modelRoot + '_' + reac.notes['abbrev']] = {'id': modelRoot + '_' + reac.notes['abbrev'], 'abbrev': reac.notes['abbrev'], 'n' : 1}
        
        reactions[reac.id] = {'id': reac.id, 'enzyme': modelRoot + '_' + reac.notes['abbrev'],  'reactants':getReacProd(reac, 'reactants'), 'products': getReacProd(reac, 'products'), 'forward_rate' : 1, 'reversed_rate': int(reac.reversibility) }
        
        
    
    return metabolites, reactions, enzymes


model = cobra.io.read_sbml_model('C:/Users/u0139894/Documents/BH_EnergyMetabolism/files/gsmms/BH_template_model.xml')

reactions = ['EX_cpd00020_e', 'rxn05938', 'rxn47229', 'Rnf', 'rxn10042']
mets,reacs,enzs = extractReactions(model, reactions, 'bh1')

with open('C:/Users/u0139894/Documents/GitHub/ChemSamp/files/full_model/enzymes.txt', 'w') as f:
    f.write('ID\tAbbrev\tn\n')
    for i in enzs:
        f.write(enzs[i]['id'] + '\t' + enzs[i]['abbrev'] + '\t' + str(enzs[i]['n']) + '\n')

with open('C:/Users/u0139894/Documents/GitHub/ChemSamp/files/full_model/metabolites.txt', 'w') as f:
    f.write('ID\tAbbrev\tn\texchange\n')
    for i in mets:
        f.write(mets[i]['id'] + '\t' + mets[i]['abbrev'] + '\t' + str(mets[i]['n']) + '\t' + str(mets[i]['exchange']) + '\n')
    
with open('C:/Users/u0139894/Documents/GitHub/ChemSamp/files/full_model/reactions.txt', 'w') as f:
    f.write('ID\tReactants\tProducts\tEnzyme\tforward_Rate\treversed_rate\n')
    for i in reacs:
        f.write(reacs[i]['id'] + '\t' + reacs[i]['reactants'] + '\t' + reacs[i]['products'] + '\t' + reacs[i]['enzyme'] + '\t' + str(reacs[i]['forward_rate']) + '\t' + str(reacs[i]['reversed_rate']) + '\n')
    
