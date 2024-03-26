#######################################################
#######################IMPORTS#########################
#######################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#######################################################
##################GLOBAL VARIABLES#####################
#######################################################

# Initialises a generator to draw values at random. (used to select random positions, initialise system and assess probabilities)
rng = np.random.default_rng()

# Create a figure and axis for the animation
fig, ax = plt.subplots()

#######################################################
######################DYNAMICS#########################
#######################################################

# Function returns the boolean if the position (i,j) has any infected neighbours
def InfectedNeighbours(state,N,i,j):

    # Indices of the 4 neighbours with state[i][j] at centre
    indices = [
        ((i-1)%N, j),   # Top
        (i, (j-1)%N),   # Left
        (i, (j+1)%N),   # Right
        ((i+1)%N, j)    # Bottom
    ]

    contact = False

    # Checks all 4 neighbours and returns True when an infected one is found
    for r, c in indices:
        if state[r][c] == 1:
            contact = True
            break

    # Returns True if been in contact with an infected individual
    return contact

# Implements the propagation rules for the SIRS Model for a full sweep of the system
def Rules(state,N,p1,p2,p3):

    for _ in range(N**2):
        # Chooses a random position in the system
        i,j = rng.choice(N,2)

        # Random float, r, is generated and used to assess against the probabilities 
        r = rng.random()

        # A Susceptible state (state[i][j] = 0) can be infected by a nearest neighbour assessed with p1
        if state[i][j] == 0:
            if InfectedNeighbours(state,N,i,j):
                if r < p1:
                    state[i][j] = 1
        # An infected state (state[i][j] = 1) can recover when assessed with p2
        if state[i][j] == 1:
            if r < p2:
                state[i][j] = 2
        # A Recovered state (state[i][j] = 2) can lose immunity and return to susceptible when assessed with p3
        if state[i][j] == 2:
            if r < p3:
                state[i][j] = 0

    # Returns updated state after a new sweep
    return state

########################################################################
###############################RUN SIM##################################
########################################################################

# Runs the Simulation with animation and option to add permanently immune individuals to the population
def RunSimA(N,sweeps,p1,p2,p3,h):

    # Initialise state of system with 1/3 probability of being any SIR compartment
    state = rng.choice([0,1,2],(N,N))

    # Initialises a sum variable that will be added to on each measurement and finally divided by the number of sweeps taken for the infection to die out or the maximum sweeps defined by the simulation params
    I_sum = 0

    # Total number of nodes in the system
    totalN = N * N

    # Randomly assigns nodes to have permanent immunity (state[i][j] == 3 has no changes in the Rules() function thus will always stay the same)
    for _ in range(int(min(h,1) * totalN)):
        i,j = rng.choice(N,2)
        state[i][j] = 3

    # Creates a multidimensional array that will store the state of the system after each sweep to be animated
    generations = np.zeros(shape=(sweeps,N,N))
    data = np.zeros(shape=(sweeps,3))

    # On first timestep save the starting state to generations array
    generations[0] = state

    # Main simulation loop runs for the specified iteration length
    for t in range(1,sweeps):

        # Performs a Sweep on the system using the SIRS model rules and returns the updated state
        state = Rules(state,N,p1,p2,p3)

        # Adds the newest state to the generations array
        print("Sweep: {0}, I: {1}".format(t,len(state[state == 1])))
        generations[t] = state
        I_sum += len(state[state == 1])

        for i in range(3):
            data[t][i] = len(state[state == i])/totalN
        
        # Breaks out of simulation loop once the infection has run its course (I = 0)
        if len(state[state == 1]) == 0:
            sweeps = t
            break

    # Uses the sum of Infected sites to calculate the average fraction. The skipped measurements to wait for equilibrium have been accounted for.
    I_mean = (I_sum/sweeps) / totalN

    return generations,I_mean,data

########################################################################

# Bootstrap Error on the Infection Variance
def Error(I_var):
    return 0

# Run Simulation to measure Infection and Variance data
def RunSimI(N,sweeps,skip,p1,p2,p3,opt):

    # Initialise state of system with 1/3 in each compartment
    state = rng.choice([0,1,2],(N,N))

    # Initialises sum variables that will be added to on each measurement to be anylysed at the end to produce the final result
    I_sum = 0
    I_sum_sq = 0

    # Total number of nodes in the system
    totalN = N * N

    # Main simulation loop runs for the specified number of sweeps
    for t in range(1,sweeps):

        # Runs a sweep over the system and updates using the rules of the SIRS model
        state = Rules(state,N,p1,p2,p3)

        # Records number of infected sites at the current sweep but only after 'skip' number of sweeps have been run
        if t > skip - 1:
            # Measures the number of Infected states at each sweep
            print("Sweep: {0:5}, p1: {1:4.3f}, p3: {2:4.3f}, I: {3}".format(t,p1,p3,len(state[state == 1])))
            I_sum += np.sum(state[state == 1])
            I_sum_sq += np.sum(state[state == 1])**2

        # Breaks out of simulation loop once the infection has run its course (I = 0)
        if len(state[state == 1]) == 0:
            sweeps = t
            break
    
    # Uses the sum of Infected sites to calculate the average fraction. The skipped measurements to wait for equilibrium have been accounted for.
    I_mean = (I_sum/(sweeps - skip - 1)) / totalN

    # Bypasses Variance and Error propagation for Infection simulations
    if opt == 'I':
        
        I_var = 0
        err_var = 0

    elif opt == 'V':
        # Uses the sum of Infected sites to calculate variance on the average fraction. The skipped measurements to wait for equilibrium have been accounted for.
        I_var = ((I_sum_sq/(sweeps - skip - 1)) - (I_sum/(sweeps - skip - 1))**2) / totalN

        # Propagates the associated error on the variance using a resampling method 
        err_var = Error(I_var)

    return I_mean,I_var,err_var

########################################################################
###############################ANIMATION################################
########################################################################

# Function called by the animation to update the frame (Each frame represents a measurement at the frequency specified at input)
def animate(t,generations,freq):
    ax.clear()
    ax.imshow(generations[t], cmap='tab20c', interpolation='nearest')
    ax.set_title(f'Sweep: {t + 1}')
    ax.set_xticks([])
    ax.set_yticks([])

########################################################################
###########################MAIN CONTROL LOOPS###########################
########################################################################

# Runs an animated plot of the simulation for the parameters inputted by the user
def mainA():
    # Initialise system parameters
    N = 100
    sweeps = 2_000

    # User inputs - outbreak fraction sets a fraction of the total nodes in the system to infected 
    p1 = float(input("p1: "))
    p2 = float(input("p2: "))
    p3 = float(input("p3: "))

    h = min(float(input("Permanent Immunity Fraction: ")),1)
    
    # Run simulation
    generations, I_mean, data = RunSimA(N,sweeps,p1,p2,p3,h)

    
    # Create and show animation
    ani = FuncAnimation(fig, animate, frames=sweeps, interval=200,fargs=(generations,1))

    # Saves the animation as a .gif file
    ani.save('Wavy.gif', writer = 'mencoder', fps=10)

    # Saves S,I,R data for each sweep to file
    f1 = open("Wave_Data.txt", 'a')
    for i in range(len(data)):
        f1.write("{0:5f}, {1:5f}, {2:5f}\n".format(data[i][0],data[i][1],data[i][2]))
    f1.close()

    plt.show()

########################################################################

# Runs the simulation for different values of p1 & p3 to find the average infected fraction associated with each parameter combination
def mainI():
    # Initialise 50x50 system, Run the simulation for 1000 sweeps wait 100 until equilibrium is reached then take a measurement of each subsequent sweep
    N = 50
    sweeps = 1_000
    skip = 50
    
    # Resolution determines the incremental factor between subsequent parameters, increment over p1 & p3 with p2 fixed at 0.5
    resolution = 0.05
    p1s = np.arange(0.3,1+resolution,resolution)
    p2 = 0.5
    p3s = np.arange(0.3,1+resolution,resolution)

    # Iterates over each combination of probability p1 & p2 and appends the Average infected sites fraction and associated p values to the data file
    for p1 in p1s:
        for p3 in p3s:
            # Run simulation
            I_mean,I_var,err_var = RunSimI(N,sweeps,skip,p1,p2,p3,'I')

            f1 = open("Infected_Data.txt", 'a')
            f1.write("{0:4.2f}, {1:4.2f}, {2}\n".format(p1,p3,abs(I_mean)))
            f1.close()

########################################################################

# Runs the simulation for different values of p1 to find the variance in infection
def mainV():
    # Initialise 50x50 system, Run the simulation for 10_000 sweeps wait 100 until equilibrium is reached then take a variance measurement at each subsequent sweep
    N = 50
    sweeps = 10_000
    skip = 100
    
    # Resolution determines the incremental factor between subsequent parameters, increment over p1 with p2 = p3 = 0.5
    resolution = 0.05
    p1s = np.arange(0.5,1+resolution,resolution)
    p2 = 0.5
    p3 = 0.5

    # Iterates over each value of probability p1 and appends the infection variance to the data file
    for p1 in p1s:
        # Run simulation
        I_mean,I_var,err_var = RunSimI(N,sweeps,skip,p1,p2,p3,'V')

        f1 = open("Variance_Data.txt", 'a')
        f1.write("{0:4.3f}, {1}, {2}\n".format(p1,I_var,err_var))
        f1.close()

########################################################################

# Runs the simulation for different values of permanent immunity to find the critical value at which infection halts propagation
def mainH():
    # Initialise 50x50 system, Run the simulation for 1000 sweeps wait 100 until equilibrium is reached then take a measurement of each subsequent sweep
    N = 50
    sweeps = 1_000
    
    # incr determines the incremental factor between subsequent fractions of permanent immunity with fixed values of p1 = p3 = 1 and p2 = 0.5
    incr = 0.005
    herds = np.arange(0,1+incr,incr)
    p1 = 1
    p2 = 0.5
    p3 = 1

    for h in herds:
        # Run simulation
        generations,I_mean = RunSimA(N,sweeps,p1,p2,p3,h)

        # Appends the Average Infected Fraction data and associated permanent immunity (h) value to the data file
        f1 = open("Herd_Immunity_Data.txt", 'a')
        f1.write("{0:4.3f}, {1}\n".format(h,abs(I_mean)))
        f1.close()

########################################################################
#############################COMMAND LINE###############################
########################################################################

# When run at the command line the user can select which simulation to run
if __name__ == "__main__":
    opt = input("Animation (A), Infected Sites (I), Variance (V), Herd Immunity (H): ").upper()
    if opt == 'I':
        mainI()
    elif opt == 'V':
        mainV()
    elif opt == 'A':
        mainA()
    elif opt == 'H':
        mainH()

########################################################################
#################################NOTES##################################
########################################################################
     
# Implement new rules:
    # S → I with probability p1 if at least one neighbour of the site to update is I; otherwise the site is unchanged.
    # I → R with probability p2.
    # R → S with probability p3.

# (i) Parameter Set Leading to an absorbing state with all sites susceptible
#     p1: 0.4
#     p2: 0.3
#     p3: 0.1
#     Or any combination of p1 & p2 with p3 = 0

# (ii) Parameter Set Leading to a dynamic equillibrium between S,I,R
#     p1: 1
#     p2: 0.5
#     p3: 1

# (iii) Parameter Set Leading to a cyclic wave of infections through the lattice
#     p1: 0.8
#     p2: 0.1
#     p3: 0.01
    