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

class PostCritical(Compartment):
    def __init__(self, name, iters):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)
    
    def Update(self, t):
        return super().Update(t)
    
    def Transmission(self):
        return super().Transmission()
    
########################################################################
################################NOTES###################################
########################################################################
    
# Link from Critical Care
# Links to Deceased and Resistant