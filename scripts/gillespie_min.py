import numpy as np
import pandas as pd

class Compound:
    '''
    Compounds are discrete entitties (molecular species) that can change their quantities based on reaction channels.
    
    To create one:
    
    c = Compound('MoleculeA', 5)
    
    After a simulation, c.n will contain the number of 'MoleculaA' that is present in the system.
    To repeat a simulation with the same object, reset the c.n to the desired starting concentration.
    
    '''
    def __init__(self, name:str, n:int):
        
        
        
        self.name = name
        self.n = n
        
        self.reactions = set()
    
    def consumeProduce(self, m:int):
        '''
        

        Parameters
        ----------
        m : int
            alter n by m. If m<0 the compound is consumed,
            otherwise it's produced. the amount of the compound is always
            greater than zero.

        Returns
        -------
        None.

        '''
        
        self.n = max(0, self.n + m)
        
        
        
    
    
            

class Reaction:
    '''
    Reaction channel defines an interaction between compounds.
    If a compound is simply consumed or degraded, provide the reactants
    or products, respectively, as empty dicts.
    '''
    def __init__(self, name:str, k:float, reactants:dict, products:dict, site=None, species=None, ecological=0):
        '''
        

        Parameters
        ----------
        name : str
            Identifier of the reaction channel
        k : float
            The rate of the reaction (the probability that 
                                      the reaction occurs in 
                                      the infinitesimal time dt)
        reactants : dict
            Reactants in the Additive interaction between compounds. 
            For example the reactants in the reaction 
            A+A--->B are represented by the dict:
                {compoundObjA:2} 
            while the reactants in the reaction A+B ---> C 
            are represented as:
                {compoundObjA:1, compoundObjB:1}
            
        products : dict
            Products in the Additive interaction between compounds. 
            For example the product in the reaction 
            A+A--->B is represented by the dict:
                {compoundObjA:1} 
            while the reactants in the reaction A+B ---> C + D 
            are represented as:
                {compoundObjC:1, compoundObjD:1}


        
        

        '''
        self.name = name
        self.k = k
        self.reactants = reactants
        self.products = products
        
        self.compounds = self.__get_compounds()
        self.update_propensity()
        self.ecological = ecological
        self.site = site
        self.species = species
    
    def pfact(self, n:int, z:int):
        '''
        Helper function to compute the propensities when the same molecule
        is consumed in a reaction channel.
        For example the reaction A+A+A--->B, the propensity is:
            nA*(nA-1)*(nA-2).k
        this function helps compute the first term of this function.
        
        

        Parameters
        ----------
        n : int
            number of molecules from a compound
        z : int
            number of times that the molecule is used
            as a reactant in the reaction channel.

        Returns
        -------
        p : int
            the first term of the propensity function.

        '''
        p=1.0
        for i in range(z):
            p*=max((n-i),0)
        return p
    
    def __get_compounds(self):
        
        compounds=set()
        
        for i in self.reactants:
            compounds.add(i)
        for i in self.products:
            compounds.add(i)
        
        
        return compounds
    
            
    def update_propensity(self):
        '''
        Update the reaction channel propensity,
        defined as the product between the quantity of reacting
        compounds and the rate. 
        (e.g. A+B-->C has the propensity defined as 
         nA*nB*k)
        When molecules of the same 
        compound interact (e.g. A+A--->B), the propensity is defined
        as nA*(nA-1)*k


        '''
        if len(self.reactants) ==0:
            
            self.propensity = self.k
            return 
        
        self.propensity = np.mean([self.pfact(i.n, int(self.reactants[i])) for i in self.reactants])*self.k
        #Applying the law of the minimum
        
    
    def react(self, mTimes):
        if self.ecological:
            a = self.site.n
            self.site.consumeProduce(-mTimes)
            b = self.site.n
            self.species.consumeProduce(a-b)
        else:
            #first check if reaction can happen
            if np.all([i.n>=mTimes*self.reactants[i] for i in self.reactants]):
                
                for i in self.reactants:
                    i.consumeProduce(-mTimes*self.reactants[i])
                    
                
                for i in self.products:
                    i.consumeProduce(mTimes*self.products[i])
        
        
                        

class ReacSystem:
    '''
    System of molecules(compounds) that interact through
    reaction channels.
    
    '''
    def __init__(self, reactions:list):
        '''
        

        Parameters
        ----------
        reactions : list
            A list of reaction objects.

        Properties
        -------
        compounds : dict
            a dict with compounds and the reactions channels
                    where they are consumed or produced.
        
        cpdIds : list 
            the names of compounds that are modelled in the system.
        
        time: np.array 
            created after a simulation is the time when the 
            system was sampled.
            
        series: np.array
            created after a simulatio. Shows the temporal evolution
            on the concetration of compounds, matched with the 
            time vector.
                    
                    

        '''
       
        self.reactions = reactions
        self.compounds = self.__get_compounds()
        self.cpdIds = [i.name for i in self.compounds]
        
    
    def __get_compounds(self):
        '''
        method to retrieve the compounds from the 
        reaction channels that are part of the system.

        Returns
        -------
        c : dict
            a dict with compounds and the reactions channels
                    where they are consumed or produced.

        '''
        
        c={}
        
        for i in self.reactions:
            
            for comp in i.compounds:
                if comp not in c:
                    c[comp]=set()
                c[comp].add(i)
                
        return c
    
    def getTau(self, tmax:float):
        '''
        Compute the time step for the direct method.

        Parameters
        ----------
        tmax : float
            The maximum time of the system. used in case the
            sum of propensities equals zero. This way we avoid
            a division by zero and no event occurs.

        Returns
        -------
        p : np.array
            probabilities associated with each reaction channel
            e.g. the normalized propensities.
        tau: float
            the time step of next reaction to be fired.

        '''
        #this loop could possibly be avoided, but it seems safe 
        #to keep it here.
        self.propensities = np.array([i.propensity for i in self.reactions])
        
        c = np.mean(self.propensities[self.propensities>0])
        
        if c==0:
            return c, tmax + 1
        
        return self.propensities/sum(self.propensities), (1/c)*np.log(1/np.random.uniform())
    
    def getTauKleaps(self, tmax:float, k:int):
        '''
        Capture the time step for the direct method. E.g. time until
        k reactions reactions occur in the system.

        Parameters
        ----------
        tmax : float
            The maximum time of the system. used in case the
            sum of propensities equals zero. This way we avoid
            a division by zero and no event occurs.
        k : int
            The number of reactions to fire per time step.

        Returns
        -------
        p: np.array
            probabilities associated with each reaction channel
            e.g. the normalized propensities used in the multinomial
            distribution.
        tau: float
            the time jump to shift the system after k reactions.

        '''
        self.propensities = np.array([i.propensity for i in self.reactions])
        c = sum(self.propensities)
        if c==0:
            return c, tmax + 1
        
        
        
        return self.propensities/c, np.random.gamma(shape = k, scale = 1/c)

    def sample(self):
        '''
        Sample the system (e.g. get the number of each compound)

        Returns
        -------
        xt : np.array
            concentration (state variable) of each compound

        '''
        return np.array([i.n for i in self.compounds])
    
    def simul(self, t:float):
        '''
        Direct method of the stochastic simulation
        algorithm.
        
        Keep in mind that in this implementation, the system
        interacts with the compounds ans reaction channels
        that are independent objects that get updated through
        this simulation. So the their concentrations change.
        If repeating a simulation, reset the concentrations
        or defined a new ReacSystem object.
        

        Parameters
        ----------
        t : float
            How long to simulate the system.


        '''
        
        
        sys_t = np.linspace(0, t, 80*60*10)
        #The system will be sampled at these times
        current_t = 0
        #a variable that keeps track of the current time
        sample_t = sys_t[0]
        # a variable the keeps track of the current sampled time
        series = np.zeros((len(sys_t), len(self.compounds)))
        # the samples
        
        
        while current_t< t:
            print(current_t)
            p, tau =self.getTau(t)
            current_t+=tau
            if current_t>t: #end simulation
                series[(sys_t>=sample_t) & (sys_t<current_t)] = self.sample()
                break
            
            if current_t>sample_t: #take a sample
                index = (sys_t>=sample_t) & (sys_t<current_t)
                series[index] = self.sample()
                sample_t = sys_t[index][-1]
                
            
            r = np.random.choice(self.reactions, p=p) #update the system
            
            r.react(1)
            for reac in self.reactions:
                reac.update_propensity()
            
            #if there are many reactions, or the system
            #is sparse, one can try to update only the
            #reactions that are affected with these ugly
            #for loops.
            # for comp in r.compounds:
            #     for reac in self.compounds[comp]:
            #         reac.update_propensity()
            
            
        
        self.time = sys_t
        self.series = series
    
    def simulKLeaps(self, t:float, k:int):
        '''
        K leaps method of the stochastic simulation
        algorithm.
        
        Perform k reactions instead of 1 per time step.
        
        Keep in mind that in this implementation, the system
        interacts with the compounds ans reaction channels
        that are independent objects that get updated through
        this simulation. So the their concentrations change.
        If repeating a simulation, reset the concentrations
        or defined a new ReacSystem object.

        Parameters
        ----------
        t : float
            How long to simulate the system.
        k : int
            How many reactions to perform per time step

        

        '''
        sys_t = np.linspace(0, t, t)
        #The system will be sampled at these times
        current_t = 0
        #a variable that keeps track of the current time
        sample_t = sys_t[0]
        # a variable the keeps track of the current sampled time
        series = np.zeros((len(sys_t), len(self.compounds)))
        # the samples
        
        
        while current_t< t:
            print(current_t)
            p, tau =self.getTauKleaps(t, k)
            current_t+=tau
            if current_t>t:#end simulation
                series[(sys_t>=sample_t) & (sys_t<current_t)] = self.sample()
                break
            
            if current_t>sample_t:#take a sample
                index = (sys_t>=sample_t) & (sys_t<current_t)
                series[index] = self.sample()
                sample_t = sys_t[index][-1]
                
            
            reactionsToFire = np.random.multinomial(k, p)#update the system
            
            
            for i,reac in enumerate(self.reactions):
                if reactionsToFire[i]>0:
                    reac.react(reactionsToFire[i])
            
            #this separate independent for loop is 
            #probably needed to assure accuracy
            for reac in self.reactions:
                reac.update_propensity()
            
            
        
        self.time = sys_t
        self.series = series



# ######Example of usage#####

# #substrate --> product

# E = Compound("enzyme", 5)
# S = Compound("substrate", 60)
# ES = Compound("enzymeSubstrate", 0)
# P = Compound("product", 0)

# r1 = Reaction("E+subs", .1, {E:1, S:1}, {ES:1})
# r2 = Reaction("E+subs diss", .1, {ES:1}, {E:1, S:1})
# r3 = Reaction("E+subs prod", .1, {ES:1}, {E:1, P:1})

# #system

# r = [r1, r2, r3]

# rsys = ReacSystem(r)

# rsys.simulKLeaps(2500,3)

# df = pd.DataFrame.from_dict({rsys.cpdIds[i]:rsys.series.T[i] for i in range(len(rsys.cpdIds))})
# df.index = rsys.time
# df.plot(title='k-leaps')


# # oscillating system
# # A = Compound("A", 10)
# # B = Compound("B", 10)

# # r1 = Reaction("1", 4*(10**-5), {A:2, B:1}, {A:3})
# # r2 = Reaction("2", 50, {}, {A:1})
# # r3 = Reaction("3", 10, {A:1}, {})
# # r4 = Reaction("4", 25, {}, {B:1})


# # #system

# # r = [r1, r2, r3, r4]

# # rsys = ReacSystem(r)

# # rsys.simul(80*60)

# # df2 = pd.DataFrame.from_dict({rsys.cpdIds[i]:rsys.series.T[i] for i in range(len(rsys.cpdIds))})
# # df2.index = rsys.time
# # df2.plot(title='k-leaps')



# 5A + B <--> C + 3D

# A = Compound("A", 600)
# B = Compound("B", 500)
# E = Compound("enzyme", 25)
# ABE = Compound("ABE", 0)
# C = Compound("C", 2)
# D = Compound("D", 5)

# r1 = Reaction("complex formation", .1, {A:5, B:1, E:1},  {ABE:1})

# r2 = Reaction("complex dissociation", .01, {ABE:1},{A:5, B:1, E:1})

# r3 = Reaction("product formation", .01, {C:1, D:3, E:1}, {ABE:1})

# r4 = Reaction("reverse", .1, {ABE:1},{C:1, D:3, E:1})
# #r2 = Reaction("product formation", 1, {ABE:1}, {A:1, B:1, E:1})
# #r3 = Reaction("product formation", 1, {ABE:1}, {C:1, D:1, E:1})
# r = [r1,r2, r3, r4]

# rsys = ReacSystem(r)

# rsys.simul(500)

# df = pd.DataFrame.from_dict({rsys.cpdIds[i]:rsys.series.T[i] for i in range(len(rsys.cpdIds))})
# df.index = rsys.time
# df.plot(title='k-leaps')

# # # 3A + B --> C + 2D
# # A = Compound("A", 600)
# # B = Compound("B", 500)
# # E = Compound("enzyme", 25)
# # ABE = Compound("complex", 5)
# # C = Compound("C", 2)
# # D = Compound("D", 2)

# # r1 = Reaction("complex formation", .00001, {A:3,B:1, E:1}, {ABE:1})
# # r2 = Reaction("dissociation", .0000001, {ABE:1}, {A:1, B:1, E:1})
# # r3 = Reaction("product formation", .1, {ABE:1}, {C:1, D:2, E:1})
# # r = [r1, r2, r3]

# # rsys = ReacSystem(r)

# # rsys.simulKLeaps(1000,10)

# # df = pd.DataFrame.from_dict({rsys.cpdIds[i]:rsys.series.T[i] for i in range(len(rsys.cpdIds))})
# # df.index = rsys.time
# # df.plot(title='k-leaps')


