########################################################################
################################IMPORTS#################################
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.integrate import simpson

########################################################################
#################################MODEL##################################
########################################################################

# Find the ppop that are moving from I to R at the current timestep using a gamma distribution
def updateR(t,newI,max_prob_R,maxT):
    # There is a maximum time a person can be infected for thus we must only search values as far back as t - maxT
    # t - maxT + 1 at t = 0 effectively runs searchT from -maxT -> 0
    searchT = np.arange(t - maxT + 1,t+1)

    # Creates an array to store the percentage of people that were infected on each day that are now recovered
    prob_R = np.zeros(len(searchT))
    i = 0
    while i < len(prob_R):
        if searchT[i] < 0:
            prob_R[i] = 0
        else:
            # Uses a gamma probability density function to find the percentage who were infected at searchT[i] who are now moving to recovered at time t
            print(searchT[i])
            print(searchT[i] - max_prob_R)
            prob_R[i] = gamma(searchT[i] - max_prob_R).pdf(searchT)[searchT[i]] * newI[searchT[i]]
            # Removes recovered proportion from the newI array
            newI[searchT[i]] = newI[searchT[i]] - prob_R[i]

        i += 1

    # Integrate over the prob_R array to determine the total amount of the population who are recovering at time t
    new_R = simpson(prob_R)
    
    return new_R 

# Find the ppop that are moving from R to S at the current timestep using a gamma distribution
def updateS(t,newR,max_prob_S):
    searchT = np.arange(0,t+1)

    # Creates an array to store the percentage of people that were resistant on each day that are now susceptible again
    prob_S = np.zeros(len(searchT))

    for i in range(len(prob_S)):
        # Uses a Gamma probability density function to find the percentage who were resistant at searchT[i] who are now moving to susceptible at time t
        prob_S[i] = gamma(searchT[i] - max_prob_S).pdf(searchT)[searchT[i]] * newR[searchT[i]]
        # Removes susceptible proportion from the newR array
        newR[searchT[i]] = newR[searchT[i]] - prob_S[i]
    
    # Integrate over the prob_S array to determine the total amount of the population who are losing immunity at time t
    new_S = simpson(prob_S)

    return new_S 

# Runs the SIR model simulation
def runModel(i_0,iters,max_prob_R,max_prob_S,maxT):
    # Initialise arrays to store the total population percentage present in that compartment at each timestep
    S = np.zeros(iters)
    I = np.zeros(iters)
    R = np.zeros(iters)

    # Initialise arrays to store the population percentage that have just entered the compartment at each timestep
    newS = np.zeros(iters)
    newI = np.zeros(iters)
    newR = np.zeros(iters)

    # Initial values normalise population to 1
    S[0] = 1 - i_0
    I[0] = i_0
    newS[0] = S[0]
    newI[0] = I[0]

    # Defines the infection rate
    beta = 0.218

    # Main loop to update values for each timestep
    for t in range(iters - 1):
        print("Day: " + str(t + 1))
        # Calculates the newly infected and recovered proportions at time t
        newI[t] = beta * S[t] * I[t]
        print(newI)
        newR[t] = updateR(t,newI,max_prob_R,maxT)
        print(newR)
        newS[t] = updateS(t,newR,max_prob_S)
        print(newS)

        # Uses the new values to update the main compartment arrays at the next timestep
        S[t+1] = S[t] + newS[t] - newI[t] 
        I[t+1] = I[t] + newI[t] - newR[t]
        R[t+1] = R[t] + newR[t] - newS[t]

    return S, I, R

########################################################################
#################################PLOT###################################
########################################################################

# Plots the Compartments values in percentage of the total population against time in days
def plot_sir_model(S, I, R, iters):
    plt.figure(figsize=(10, 6))
    plt.plot(range(iters), 100*S, label='Susceptible')
    plt.plot(range(iters), 100*I, label='Infected')
    plt.plot(range(iters), 100*R, label='Resistant')
    plt.xlabel('Time (Days)')
    plt.ylabel('% of Population')
    plt.title('SIRS Model Simulation')
    plt.legend()
    plt.show()

########################################################################
#############################COMMAND LINE###############################
########################################################################
    
# Allows code to be run from the command line, it takes input parameters then runs and plots the results of the simulation
if __name__ == "__main__":
    # Input initialisation parameters: i_0 - initial % of population infected
    i_0 = float(input("Initial Percentage of population infected [0,100]: ")) / 100
    iters = int(input("Number of Timesteps: ")) + 1
    max_prob_R = int(input("How many days after initial infection do most people recover? ")) + 1
    maxT = int(input("Maximum time someone can be infected for? "))
    max_prob_S = int(input("How many days after recovery do most people lose their immunity? ")) + 1

    # Run SIR model
    S, I, R = runModel(i_0,iters,max_prob_R,max_prob_S,maxT)

    # Check population is conserved to an acceptible tolerance
    if (S[-1]+I[-1]+R[-1]).round(6) == 1:
        print("Population is conserved")
    
    # Plot data
    plot_sir_model(S, I, R, iters)

########################################################################
################################NOTES###################################
########################################################################
    
