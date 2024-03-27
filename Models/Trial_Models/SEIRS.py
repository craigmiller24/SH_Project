import numpy as np
import matplotlib.pyplot as plt

def seirs_model(beta, gamma, alpha, mu, omega, sigma, days):
    S = np.zeros(days)
    E = np.zeros(days)
    I = np.zeros(days)
    R = np.zeros(days)

    # Initial values normalise population to 1
    S[0] = 0.999
    E[0] = 0.001
    I[0] = 0
    R[0] = 0

    # SIIR model simulation
    for t in range(1,days):
        S[t] = max(0,S[t-1] + (mu - (beta * I[t-1] * S[t-1]) + (omega * R[t-1]) - (mu * S[t-1])))
        E[t] = max(0,E[t-1] + ((beta * I[t-1] * S[t-1]) - (sigma * E[t-1]) - (mu * E[t-1])))
        I[t] = max(0,I[t-1] + ((sigma * E[t-1]) - (gamma * I[t-1]) - ((mu + alpha) * I[t-1])))
        R[t] = max(0,R[t-1] + ((gamma * I[t-1]) - (omega * R[t-1]) - (mu * R[t-1])))

    print(100*E[-1])
    print(100*I[-1])

    return S,E,I,R

def plot_seirs_model(susceptible, exposed, infected, recovered, days):
    plt.figure(figsize=(10, 6))
    plt.yscale('log')
    #plt.plot(range(days), 100*susceptible, label='Susceptible')
    plt.plot(range(days), 100*exposed, label='Exposed')
    plt.plot(range(days), 100*infected, label='Infected')
    #plt.plot(range(days), 100*recovered, label='Recovered')
    plt.xlabel('Time (Days)')
    plt.ylabel("ln(% of Population)")
    plt.title('SEIRS Model Simulation')
    plt.legend()
    plt.show()

# Parameters
beta = float(input("Infection Rate: "))             # Infection rate 0.218 => Estimated R_0 = 2.4
gamma = 1/float(input("Recovery Rate: "))           # Recovery rate 11 days
alpha = float(input("I to R Rate: "))               # 0 deaths due to infection
mu = 1/float(input("Birth/Death Rate: "))           # 29_528 days (80.9 years - UK life expectancy 2023)
omega = 1/float(input("Immunity Loss Rate: "))      # 365 (1 year)
sigma = 1/float(input("Latency Rate: "))            # 5.5 days
days = int(input("Simulation Period: "))            # Simulated time period

# Run the SEIR model
susceptible, exposed, infected, recovered = seirs_model(beta, gamma, alpha, mu, omega, sigma, days)

# Plot the results
plot_seirs_model(susceptible, exposed, infected, recovered, days)
