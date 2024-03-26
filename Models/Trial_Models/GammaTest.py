import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.stats import norm
from scipy.integrate import simpson

# Gamma Dist(a,x) = x**(a-1) * exp(-x)
# a - 1 is the position of the peak of the function
t = np.arange(0,100,1)

a = float(input("max: "))
pgamma = gamma(a).pdf(t)
pnorm = norm.pdf(t, loc= a, scale=0.15*a)

fig = plt.figure()

plt.plot(t,pgamma,label="Gamma: a = " + str(a))
plt.plot(t,pnorm,label="Gaussian: a = " + str(a))
plt.xlabel("Time (Days)")
plt.ylabel("Probability of Transition, P(t)")
plt.title("Probability Distributions")
plt.legend()

plt.show()