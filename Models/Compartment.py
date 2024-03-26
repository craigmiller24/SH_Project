########################################################################
################################IMPORTS#################################
########################################################################

# Modules
import numpy as np

########################################################################
#################################CLASS##################################
########################################################################

# Parent class for all compartments in the model
class Compartment(object):
    # Initialise fields
    def __init__(self,name,iters):
        self.name = name
        self.vals = np.zeros(iters)
        self.new_val = np.zeros(iters)

    # Called at each iteration of the simulation to find out the proportion of the population that are entering this compartment at that timestep
    def Update(self,t):
        return self.new_val[t]
    
    # Called from the update function to determine how the people leaving the compartment are to be distributed amongst the other compartments
    # Need a field that contains all of a classes connections and their associated probability
    #   Example: Infected people can move to one of resistant, hospitalised or deceased this could have a probability ratio of 90:9:1
    def Transmission(self):
        return None
    
########################################################################
################################NOTES###################################
########################################################################