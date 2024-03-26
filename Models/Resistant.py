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

class Resistant(Compartment):
    def __init__(self,name,iters,max_prob,maxT):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)
        self.max_prob = max_prob
        self.maxT = maxT

    # Find the proportion of the population that are moving from I to R at the current timestep using a gamma distribution
    def Update(self,t,newI):
        # There is a maximum time a person can be infected for thus we must only search values as far back as t - maxT
        # t - maxT + 1 at t = 0 effectively runs searchT from -maxT -> 0
        searchT = np.arange(t - self.maxT + 1,t+1)

        # Creates an array to store the percentage of people that were infected on each day that are now recovered
        prob_R = np.zeros(len(searchT))
        i = 0
        while i < len(prob_R):
            if searchT[i] < 0:
                prob_R[i] = 0
            else:
                # Uses a gamma probability density function to find the percentage who were infected at searchT[i] who are now moving to recovered at time t
                prob_R[i] = gamma(self.max_prob).pdf(searchT[i]) * newI[searchT[i]]
                # Removes recovered proportion from the newI array
                newI[searchT[i]] = newI[searchT[i]] - prob_R[i]

            i += 1

        # Integrate over the prob_R array to determine the total amount of the population who are recovering at time t
        self.new_val[t] = simpson(prob_R)

        return self.new_val[t] 
    
    def Transmission(self):
        return None
    
########################################################################
################################NOTES###################################
########################################################################
    
# Links from All compartments
# Links to Deceased and Susceptible