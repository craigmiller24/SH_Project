########################################################################
################################IMPORTS#################################
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma

########################################################################
#################################MODEL##################################
########################################################################

# Find the ppop that are moving from I to R at the current timestep using a gamma distribution
def updateR(currT,newI,max_prob_R):

    # Only need to consider 95% of the gamma distributions
    maxT_Infected = int(round(gamma.ppf(0.95, max_prob_R)))

    # There is a maximum time a person can be infected for, thus we must only search values as far back as 95% confidence
    searchbackT = np.arange(currT - maxT_Infected + 1,currT + 1)

    # Creates an array to store the percentage of people that were infected on each day that are now recovered
    prob_R = np.zeros(maxT_Infected)
    day = 0

    while day < maxT_Infected:
        if currT < maxT_Infected:
            prob_R[day] = 0
        else:
            # Uses a gamma probability density function to find the percentage who were infected at searchT[i] who are now moving to recovered at time t
            prob_R[day] = gamma(currT - max_prob_R).pdf(searchbackT)[day] * newI[day]
            print(gamma(currT - max_prob_R).pdf(searchbackT)[day]*newI[day])
            print(prob_R)
            # Removes recovered proportion from the newI array
            newI[day] = newI[day] - prob_R[day]

            # On last day of search all individuals still infected must be moved to recovered to conserve population
            if day == maxT_Infected - 1:
                newI[day] -= newI[currT-maxT_Infected]
                newI[day] = min(1,max(newI[day],0))
                newI[currT-maxT_Infected] = 0
                prob_R[day] += newI[day]

        day += 1

    # Integrate over the prob_R array to determine the total amount of the population who are recovering at time t
    recovered = np.sum(prob_R)

    return recovered 

# Find the ppop that are moving from R to S at the current timestep using a gamma distribution
def updateS(currT,newR,max_prob_S):
    # Only need to consider 95% of the gamma distributions
    maxT_Immunity = int(round(gamma.ppf(0.95, max_prob_S)))

    # There is a maximum time a person can be immune for, thus we must only search values as far back as 95% confidence
    searchbackT = np.arange(currT - maxT_Immunity + 1,currT + 1)

    # Creates an array to store the percentage of people that were resistant on each day that are now susceptible again
    prob_S = np.zeros(maxT_Immunity)
    day = 0

    while day < maxT_Immunity:
        if currT < maxT_Immunity:
            prob_S[day] = 0
        else:
            # Uses a Gamma probability density function to find the percentage who were resistant at searchT[i] who are now moving to susceptible at time t
            prob_S[day] = gamma(currT - max_prob_S).pdf(searchbackT)[day] * newR[day]

            # Removes newly susceptible proportion from the newR array
            newR[day] = newR[day] - prob_S[day]
            
            # On last day of search all individuals still with immunity must be moved to susceptible to conserve population
            if day == maxT_Immunity - 1:
                newR[day] -= newR[currT-maxT_Immunity]
                newR[day] = min(1,max(newR[day],0))
                newR[currT-maxT_Immunity] = 0
                prob_S[day] += newR[day]

        day += 1
    
    # Integrate over the prob_S array to determine the total amount of the population who are losing immunity at time t
    susceptible = np.sum(prob_S)

    return susceptible 

# Runs the SIR model simulation
def runModel(i_0,iters,max_prob_R,max_prob_S,R0):
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
    beta = R0 / max_prob_R

    # Main loop to update values for each timestep
    for t in range(iters - 1):
        # Calculates the newly infected and recovered proportions at time t
        newI[t] = beta * S[t] * I[t]
        newR[t] = updateR(t,newI,max_prob_R)
        newS[t] = updateS(t,newR,max_prob_S)

        # Uses the new values to update the main compartment arrays at the next timestep
        ds = newS[t] - newI[t]
        di = newI[t] - newR[t]
        dr = newR[t] - newS[t]

        #print(ds)
        #print(di)
        #print(dr)

        S[t+1] = S[t] + ds
        I[t+1] = I[t] + di
        R[t+1] = R[t] + dr

    return S, I, R

########################################################################
#################################PLOT###################################
########################################################################

# Plots the Compartments values in percentage of the total population against time in days
def plot_sir_model(S, I, R, iters):
    plt.figure(figsize=(12, 6))
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
    R0 = float(input("Basic Reproduction Number, R_0: "))
    max_prob_R = int(input("How many days after initial infection do most people recover? ")) + 1
    max_prob_S = int(input("How many days after recovery do most people lose their immunity? ")) + 1
    
    # Run SIR model
    S, I, R = runModel(i_0,iters,max_prob_R,max_prob_S,R0)

    # Check population is conserved to an acceptible tolerance
    if (S[-1]+I[-1]+R[-1]).round(6) == 1:
        print("Population is conserved")
    
    # Plot data
    plot_sir_model(S, I, R, iters)

########################################################################
################################NOTES###################################
########################################################################
    
