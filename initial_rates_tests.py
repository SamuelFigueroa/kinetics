#!/usr/bin/env python
# coding: utf-8

# ## Determining the rate law of chemical reactions
# ## Tests

# In[1]:


from initial_rates import *
from matplotlib import pyplot as plt


#  #### Example 1:
#  Consider the concentration of $[\mathrm{N_{2}O_{5}}]$ as a function of time for the following reaction at 318 K.
#  <center>$2\space\mathrm{N_{2}O_{5}(g)} \longrightarrow 2\space\mathrm{NO_{2}(g)} + \frac{1}{2}\space\mathrm{O_{2}(g)}$</center>
# <center>
#     $\begin{array} {r|r}
#     t/\textrm{min} & [\mathrm{N_{2}O_{5}}]/10^{-2}M
#     \\\hline
#     0 & 1.24 \\
#     10 & 0.92 \\
#     20 & 0.68 \\
#     30 & 0.50 \\
#     40 & 0.37 \\
#     50 & 0.28 \\
#     60 & 0.20 \\
#     70 & 0.15 \\
#     80 & 0.11 \\
#     90 & 0.08 \\
#     100 & 0.06 \\
#     \end{array}
# $</center>

# In[2]:


# Test compInitialRate
stoichiometric_coefficient = 2

# Create array of time measurements in units of seconds.
time_points = numpy.multiply(60,[*range(0,110,10)])

# Create array of concentrations in units of moles per liter.
concentrations = numpy.multiply(1e-2, [1.24, 0.92, 0.68, 0.50, 0.37, 0.28, 0.20, 0.15, 0.11, 0.08, 0.06])

# Calculate initial rate.
initial_rate = compInitialRate(stoichiometric_coefficient, time_points, concentrations, reaction_percentage = 10.0)

# Generate initial rate curve.
species_is_reactant = True
initial_rate_curve = concentrations[0] + numpy.multiply(
    (-1 if species_is_reactant else 1) * stoichiometric_coefficient * initial_rate,
    time_points)

# Plot kinetic data and the initial rate curve.
plt.title("Initial rate of reaction from kinetic data") 
plt.xlabel("t/s") 
plt.ylabel("[N2O5]/M") 
plt.plot(time_points, concentrations, "-b", label = "[N2O5] at t")
plt.plot(time_points, initial_rate_curve,  "-r", label = "Initial rate")
plt.legend(loc = "upper right")
plt.ylim(0, concentrations[(0 if species_is_reactant else -1)])
plt.xlim(0, time_points[-1])
plt.show()
print("initial_rate = ", initial_rate, "M/s")


# #### Example 2:
# Consider the following initial rate data for the reaction.
# <center>$2\space\mathrm{NO_{2}(g)}+ \mathrm{F_{2}(g)} \longrightarrow 2\space\mathrm{NO_{2}F(g)}$</center><br/>
# <center>
#     $\begin{array} {r|r|r|r}
#     \textrm{Run} & [\mathrm{NO_{2}}]_{0}/M & [\mathrm{F_{2}}]_{0}/M & v(0)/M\cdot s^{-1}
#     \\\hline
#     1 & 1.15 & 1.15 & 6.12\space\mathrm{x}\space10^{-4} \\
#     2 & 1.72 & 1.15 & 1.36\space\mathrm{x}\space10^{-3} \\
#     3 & 1.15 & 2.30 & 1.22\space\mathrm{x}\space10^{-3} \\
#     \end{array}
# $</center>

# In[3]:


# Test compOrder
initial_concs = [[1.15, 1.15],
                 [1.72, 1.15],
                 [1.15, 2.3]]
initial_rates = [6.12e-4, 1.36e-3, 1.22e-3]
m_NO2 = compOrder(initial_rates[0], initial_rates[1], initial_concs[0][0], initial_concs[1][0])
m_F2 = compOrder(initial_rates[0], initial_rates[2], initial_concs[0][1], initial_concs[2][1])

# Check that order of NO2 is 2
assert abs(m_NO2-2.0) < 0.05
# Check that order of F2 is 1
assert abs(m_F2-1.0) < 0.05


# In[4]:


# Test compRateConstant
orders = [m_NO2, m_F2]
k = compRateConstant(initial_concs, initial_rates, orders)

# Check that k is 4.04e-4
assert (abs(k-4.04e-4)) < 1e-6


# In[5]:


# Test compLeastSquaresRateLaw
[k, m_NO2, m_F2] = compLeastSquaresRateLaw(initial_concs, initial_rates)

# Check that order of NO2 is 2
assert abs(m_NO2-2.0) < 0.05
# Check that order of F2 is 1
assert abs(m_F2-1.0) < 0.05
# Check that k is 4.04e-4
assert (abs(k-4.04e-4)) < 1e-6

