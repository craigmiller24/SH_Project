########################################################################
################################IMPORTS#################################
########################################################################

# Modules
import numpy as np
from scipy.stats import gamma
from scipy.integrate import simpson

# Classes
from Compartment import Compartment

########################################################################
#################################CLASS##################################
########################################################################

# max_capacity defines the critical care capacity which determines whether lockdown measures are required which then triggers a reduction in beta
class CriticalCare(Compartment):
    def __init__(self, name, iters, max_capacity):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)
        self.max_capacity = max_capacity
    
    def Update(self, t):
        return super().Update(t)
    
    def Transmission(self):
        return super().Transmission()
    
########################################################################
################################NOTES###################################
########################################################################
    
# Link from Hospitalised
# Links to Deceased, Post Critical and Resistant