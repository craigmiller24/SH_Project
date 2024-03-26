import numpy as np
import matplotlib.pyplot as plt

def sir_model(beta, gamma, i_0, iterations):
    # Arrays to store the results over time
    S = np.zeros(iterations + 1)
    I = np.zeros(iterations + 1)
    R = np.zeros(iterations + 1)
    
    # Initial values normalise population to 1
    S[0] = 1 - i_0
    I[0] = i_0

    # SIRS model simulation
    for i in range(iterations):
        new_I = beta * S[i] * I[i]
        new_R = gamma * I[i]

        S[i+1] = S[i] - new_I 
        I[i+1] = I[i] + new_I - new_R
        R[i+1] = R[i] + new_R

    return S, I, R

def plot_sir_model(S, I, R, iterations):
    plt.figure(figsize=(10, 6))
    plt.plot(range(iterations + 1), S, label='Susceptible')
    plt.plot(range(iterations + 1), I, label='Infected')
    plt.plot(range(iterations + 1), R, label='Recovered')
    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.title('SIR Model Simulation')
    plt.legend()
    plt.show()

# Parameters
beta = float(input("Infection Rate: "))             # Infection rate
gamma = float(input("Recovery Rate: "))             # Recovery rate
i_0 = float(input("Initial outbreak fraction: "))   # Initial outbreak fraction
timestep = float(input("Simulation Timestep"))      # Simulation timestep
days = int(input("Simulation Period: "))            # Simulated time period
iterations = round(days * 24/timestep)

# Run the SIR model
S, I, R = sir_model(beta, gamma, i_0, iterations)

# Plot the results
plot_sir_model(S, I, R, iterations)