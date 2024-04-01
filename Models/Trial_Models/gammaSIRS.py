import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma

# Define parameters
population = 10000   # population
beta = 0.1          # transmission rate
gamma_shape_R = 7   # shape parameter of gamma distribution for recovery rate
gamma_shape_S = 80  # shape parameter of gamma distribution for immunity rate

# Define initial conditions
I0 = 1
R0 = 0
S0 = population - I0

# Time vector
t = np.arange(0, 365)

# Define SIRS model equations
def SIRS(y, t, beta, gamma_shape_R, gamma_shape_S):
    S, I, R = y
    dSdt = gamma.pdf(t, gamma_shape_S,loc=0) * R - beta * S * I
    dIdt = beta * S * I - gamma.pdf(t, gamma_shape_R,loc=0) * I
    dRdt = gamma.pdf(t, gamma_shape_R,loc=0) * I - gamma.pdf(t, gamma_shape_S,loc=0) * R
    
    return dSdt, dIdt, dRdt


# Solve the SIRS equations
from scipy.integrate import odeint
y0 = S0, I0, R0
args = (beta, gamma_shape_R, gamma_shape_S)
solution = odeint(SIRS, y0, t, args=args)


# Plot the results
plt.figure(figsize=(15,5))
plt.plot(t, 100 * solution[:,0]/population, label='Susceptible')
plt.plot(t, 100 * solution[:,1]/population, label='Infected')
plt.plot(t, 100 * solution[:,2]/population, label='Resistant')
plt.xlabel('Time (Days)')
plt.ylabel('% of Population')
plt.title('SIRS Model (with Gamma)')
plt.legend()
plt.show()
