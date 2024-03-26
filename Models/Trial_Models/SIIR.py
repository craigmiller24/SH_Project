import numpy as np
import matplotlib.pyplot as plt

def siir_model(beta, gamma, alpha, i_0, days):
    # Initial values normalise population to 1
    susceptible = 1 - i_0
    infected_1 = i_0
    infected_2 = 0
    recovered = 0

    # Lists to store the results over time
    susceptible_list = [susceptible]
    infected_list = [infected_1 + infected_2]
    recovered_list = [recovered]

    # SIIR model simulation
    for day in range(days):
        new_infected_s1 = beta * susceptible * (infected_1 + infected_2)
        new_infected_s2 = alpha * infected_1
        new_recovered = gamma * infected_2

        susceptible -= new_infected_s1 
        infected_1 += new_infected_s1 - new_infected_s2
        infected_2 += new_infected_s2 - new_recovered
        recovered += new_recovered

        # Ensure values don't go below 0
        susceptible = max(susceptible, 0)
        infected_1 = max(infected_1, 0)
        infected_2 = max(infected_2, 0)
        recovered = max(recovered, 0)

        susceptible_list.append(susceptible)
        infected_list.append(infected_1 + infected_2)
        recovered_list.append(recovered)
        print(susceptible+infected_1+infected_2+recovered)
    return susceptible_list, infected_list, recovered_list

def plot_sirs_model(susceptible, infected, recovered, days):
    plt.figure(figsize=(10, 6))
    plt.plot(range(days + 1), susceptible, label='Susceptible')
    plt.plot(range(days + 1), infected, label='Infected')
    plt.plot(range(days + 1), recovered, label='Recovered')
    plt.xlabel('Days')
    plt.ylabel('Population')
    plt.title('SIIR Model Simulation')
    plt.legend()
    plt.show()

# Parameters
beta = float(input("Infection Rate: "))             # Infection rate
gamma = float(input("Recovery Rate: "))             # Recovery rate
alpha = float(input("Latency Rate: "))              # Rate at which recovered individuals become susceptible again
i_0 = float(input("Initial outbreak fraction: "))   # Initial outbreak fraction
days = int(input("Simulation Period: "))            # Simulated time period

# Run the SIIR model
susceptible, infected, recovered = siir_model(beta, gamma, alpha, i_0, days)

# Plot the results
plot_sirs_model(susceptible, infected, recovered, days)
