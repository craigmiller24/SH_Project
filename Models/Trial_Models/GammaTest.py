import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma
from scipy.stats import norm
from scipy.integrate import simpson

# Gamma Dist(a,x) = x**(a-1) * exp(-x)
# a - 1 is the position of the peak of the function
maxprob = 70
currT = 365
maxT = gamma.ppf(0.95, currT - maxprob + 1)
print(maxT)
searchbackT = np.arange(currT - maxT + 1,currT + 1)

pgamma = gamma(currT - maxprob + 1).pdf(searchbackT)

fig = plt.figure()

plt.plot(searchbackT,pgamma,label="Gamma: a = " + str(maxprob))
plt.xlabel("Time (Days)")
plt.ylabel("Probability of Transition, P(t)")
plt.title("Probability Distributions")
plt.legend()

plt.show()