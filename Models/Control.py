########################################################################
################################IMPORTS#################################
########################################################################

# Modules
import numpy as np
import matplotlib.pyplot as plt

# Classes
from Susceptible import Susceptible
from Infected import Infected
from Resistant import Resistant
from Hospitalised import Hospitalised
from CriticalCare import CriticalCare
from PostCritical import PostCritical
from Deceased import Deceased

########################################################################
#################################MODEL##################################
########################################################################

# Runs the SIR model simulation
def RunModel(i_0,iters,S,I,R,H,CC,PC,D):

    # Initial values normalise population to 1
    S.vals[0] = 1 - i_0
    I.vals[0] = i_0
    S.new_val[0] = S.vals[0]
    I.new_val[0] = I.vals[0]

    # Main loop to update values for each timestep
    for t in range(iters-1):
        # Calculates the newly infected and recovered proportions at time t
        I.new_val[t] = I.Update(t,S.vals)
        H.new_val[t] = R.Update(t,I.new_val)
        CC.new_val[t] = R.Update(t,H.new_val)
        PC.new_val[t] = R.Update(t,CC.new_val)
        
        new_rs = (I.new_val[t] + H.new_val[t] + CC.new_val[t] + PC.new_val[t])
        R.new_val[t] = R.Update(t,I.new_val)

        S.new_val[t] = S.Update(t,R.new_val)
        
        new_ds = mu * (S.new_val[t] + I.new_val[t] + H.new_val[t] + CC.new_val[t] + PC.new_val[t])
        D.new_val[t] = D.Update(t,new_ds)

        #print('Day {0:1}: dS = {1:1.3e}, dI = {2:1.3e}, dR = {3:1.3e}'.format(t, S.new_val[t]-I.new_val[t], I.new_val[t]-R.new_val[t], R.new_val[t]-S.new_val[t]))
        
        # Uses the new values to update the main compartment arrays at the next timestep
        S.vals[t+1] = S.vals[t] + S.new_val[t] - I.new_val[t] 
        I.vals[t+1] = I.vals[t] + I.new_val[t] - R.new_val[t]
        R.vals[t+1] = R.vals[t] + R.new_val[t] - S.new_val[t]

    return S, I, R

########################################################################
#################################PLOT###################################
########################################################################

# Plots the Compartments values in percentage of the total population against time in days
def plot_sirs_model(compartments,iters,cc):
    plt.figure(figsize=(10, 6))

    # Plots the data curves for each compartment
    for c in compartments:
        plt.plot(np.arange(iters), c.vals, label=c.name)

    # Plots a red dashed line at the level of the maximum ICU capacity
    plt.plot(np.arange(iters),[cc.max_capacity]*iters,'r--',label='Critical Care Capactity')

    plt.xlim(0,iters)
    plt.ylim(0,1)
    plt.xlabel('Time (Days)')
    plt.ylabel('% of Population')
    plt.title('SIR Model Simulation')
    plt.legend()
    plt.show()

########################################################################
#############################COMMAND LINE###############################
########################################################################
    
# Allows code to be run from the command line, it takes input parameters then runs and plots the results of the simulation
if __name__ == "__main__":
    # Input initialisation parameters: i_0 - initial % of population infected
    i_0 = float(input("Initial Percentage of population infected [0,100]: ")) / 100
    beta = float(input("Initial Infection Rate: "))
    iters = int(input("Number of Timesteps: ")) + 1
    max_prob_R = int(input("How many days after initial infection do most people recover? ")) + 1
    maxT = int(input("Maximum time someone can be infected for? "))
    max_prob_S = int(input("How many days after recovery do most people lose their immunity? ")) + 1
    max_capacity = float(input("Maximum Critical Care Capacity: "))

    # Use inputs to create class instances of each compartment
    S = Susceptible('Susceptible',iters,max_prob_S)
    I = Infected('Infected',iters,beta)
    R = Resistant('Resistant',iters,max_prob_R,maxT)
    H = Hospitalised('Hospitalised',iters)
    CC = CriticalCare('Critical Care',iters,max_capacity)
    PC = PostCritical('Post Critical',iters)
    D = Deceased('Deceased',iters)

    # Runs model and returns a list of the compartments
    compartments = RunModel(i_0,iters,S,I,R,H,CC,PC,D)
    
    # Plot data
    plot_sirs_model(compartments,iters,CC)

########################################################################
################################NOTES###################################
########################################################################
    
# This file is the main control file of the model, it takes inputs and initialises all of the compartments as class instances and runs the model using their attributes. It then plots the results
# Set BCs s.t. S, I & R stay within the range: [0,1] but conserve the population (e.g. sum(pop) == 1)
# Model needs generalised beyond SIR s.t all compartments are in a compartment array and the model is ran from that and thus no need to explicitly name specitic objects