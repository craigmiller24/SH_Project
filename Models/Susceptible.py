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

class Susceptible(Compartment):
    # Initialises variables
    def __init__(self,name,iters,max_prob):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)
        self.max_prob = max_prob

    # Find the proportion of the population that are moving from R to S at the current timestep using a gamma distribution
    def Update(self,t,newR):
        searchT = np.arange(0,t+1)

        # Creates an array to store the percentage of people that were resistant on each day that are now susceptible again
        prob_S = np.zeros(len(searchT))

        for i in range(len(prob_S)):
            # Uses a Gamma probability density function to find the percentage who were resistant at searchT[i] who are now moving to susceptible at time t
            prob_S[i] = gamma(self.max_prob).pdf(searchT[i]) * newR[searchT[i]]
            # Removes susceptible proportion from the newR array
            newR[searchT[i]] = newR[searchT[i]] - prob_S[i]
        
        # Integrate over the prob_S array to determine the total amount of the population who are losing immunity at time t
        self.new_val[t] = simpson(prob_S)

        return self.new_val[t] 
    
    def Transmission(self):
        return None
    
########################################################################
################################NOTES###################################
########################################################################
    
# Link from Resistant
# Links to Deceased, Infected and Resistant (via Vaccine)