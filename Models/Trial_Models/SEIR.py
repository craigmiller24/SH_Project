import numpy as np
import matplotlib.pyplot as plt

def seir_model(beta, gamma, alpha, i_0, days):
    # Initial values normalise population to 1
    susceptible = 1 - i_0
    exposed = i_0
    infected = 0
    recovered = 0

    # Lists to store the results over time
    susceptible_list = [susceptible]
    infected_list = [exposed + infected]
    recovered_list = [recovered]

    # SIIR model simulation
    for day in range(days):
        new_exposed = beta * susceptible * (exposed + infected)
        new_infected = alpha * exposed
        new_recovered = gamma * infected

        susceptible -= new_exposed 
        exposed += new_exposed - new_infected
        infected += new_infected - new_recovered
        recovered += new_recovered

        # Ensure values don't go below 0
        susceptible = max(susceptible, 0)
        exposed = max(exposed, 0)
        infected = max(infected, 0)
        recovered = max(recovered, 0)

        susceptible_list.append(susceptible)
        infected_list.append(exposed + infected)
        recovered_list.append(recovered)

    return susceptible_list, infected_list, recovered_list

def plot_seir_model(susceptible, infected, recovered, days):
    plt.figure(figsize=(10, 6))
    plt.plot(range(days + 1), susceptible, label='Susceptible')
    plt.plot(range(days + 1), infected, label='Infected')
    plt.plot(range(days + 1), recovered, label='Recovered')
    plt.xlabel('Time (Days)')
    plt.ylabel('Population')
    plt.title('SEIR Model Simulation')
    plt.legend()
    plt.show()

# Parameters
beta = float(input("Infection Rate: "))             # Infection rate
gamma = float(input("Recovery Rate: "))             # Recovery rate
alpha = float(input("Latency Rate: "))              # Rate at which recovered individuals become susceptible again
i_0 = float(input("Initial outbreak fraction: "))   # Initial outbreak fraction
days = int(input("Simulation Period: "))            # Simulated time period

# Run the SEIR model
susceptible, infected, recovered = seir_model(beta, gamma, alpha, i_0, days)

# Plot the results
plot_seir_model(susceptible, infected, recovered, days)
