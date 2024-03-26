########################################################################
################################IMPORTS#################################
########################################################################

# Modules
import numpy as np

# Classes
from Compartment import Compartment

########################################################################
#################################CLASS##################################
########################################################################

class Infected(Compartment):
    # Initialises variables
    def __init__(self,name,iters,beta):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)
        self.beta = beta

    # Called at each iteration of the simulation to find out the proportion of the population that are entering this compartment at that timestep
    def Update(self,t,S):
        self.new_val[t] = self.beta * S[t] * self.vals[t]
        return self.new_val[t]
    
    def Transmission(self):
        return None
    
########################################################################
################################NOTES###################################
########################################################################
    
# Link from Susceptible
# Links to Deceased, Hospitalised and Resistant