#!/usr/bin/env python
# coding: utf-8

# # Determining the rate law of chemical reactions
# 
# #### I. Definition of a rate law
# For a general chemical reaction described by
# 
# <center>$\nu_{A}A + \nu_{B}B \longrightarrow \nu_{Y}Y + \nu_{Z}Z,$</center><br/>
# 
# the *rate of reaction*, $v(t)$, is defined as
# 
# <center> $v(t) = -\frac{1}{\nu_{A}}\frac{d[A]}{dt} = -\frac{1}{\nu_{B}}\frac{d[B]}{dt} = \frac{1}{\nu_{Y}}\frac{d[Y]}{dt} = \frac{1}{\nu_{Z}}\frac{d[Z]}{dt} = \frac{1}{V}\frac{d\xi}{dt},$ </center><br/>
# 
# where $V$ is the volume of the system and $\xi$ is the *extent of the reaction*. The extent of the reaction is defined such that
# 
# <center>$n_{A}(t) = n_{A}(0) - \nu_{A}\xi(t),$</center><br/>
# 
# where $n_{A}$ is the number of moles of reactant $A$.
# 
# In general, $v(t)$ is related to the concentrations of the various chemical species present at time $t$.
# 
# The differential equations that describe this relationship are known as *rate laws*.
# 
# In the context of the above general reaction, a rate law may be of the form
# 
# <center>$v(t) = -\frac{1}{\nu_{A}}\frac{d[A]}{dt} = k[A]^{m_{A}}[B]^{m_{B}},$</center><br/>
# 
# where the exponents, known as *orders*, $m_{A}$ and $m_{B}$ are constants. The proportionality constant $k$ is called the *rate constant*.
# 
# #### II. Key assumptions and observations
# 1. $v(t)$ is a positive quantity.
# 2. In the expressions used above to define the reaction rate, it is assumed that the volume of the system $V$ is constant.
# 3. The order of a reactant in a rate law often differs from its stoichiometric coefficient in the balanced chemical reaction equation.
# 4. The units of the rate constant generally vary among different rate laws. This ensures that $v(t)$ has units of $\frac{concentration}{time}$.
# 5. Many rate laws are not of the form given above. Chemical reactions that occur as multi-step processes may have more complicated rate laws. An example of these may be of the form
# 
# <center> $v(t) = -\frac{1}{\nu_{A}}\frac{d[A]}{dt} = \frac{k'[A]^{m_{A}}[B]^{m'_{B}}}{1+k''[Y]^{m_{Y}}[B]^{m''_{B}}}$. </center><br/>
# 

# ## Experimental Methods
# #### I. Method of initial rates
# 
# It is not possible to directly measure the instantaneous (differential) rate at a moment $t'$, $\frac{d[A]}{dt}\Bigr|_{\substack{t=t'}}$. We can measure the change in concentration over a finite period of time $\Delta t$ and use it to approximate the instantaneous rate at any particular moment $t'$ within the time period.
# 
# 
# In most cases, the concentration of the species $A$, $[A]$, is measured a certain number of times from $t=0$ to $t=t'$. From these measurements, the instantaneous rate $\frac{d[A]}{dt}\Bigr|_{\substack{t=0}}$ can be estimated and then related to $v(0)$. 
# 
# <center>$v(0) = -\frac{1}{\nu_{A}}\frac{d[A]}{dt}\Bigr|_{\substack{t=0}} = k[A]_{0}^{m_{A}}[B]_{0}^{m_{B}}$</center><br/>
# 
# The *initial* rate, $v(0)$, is studied because the concentrations of each of the chemical species are known at $t=0$. At that instant, $[A]$ is $[A]_{0}$ and $[B]$ is $[B]_{0}$. If no product is present at $t=0$, a reverse reaction will not affect the measurement.
# 
# Consider then two experiments carried out to determine $k$, $m_{A}$ and $m_{B}$. The initial concentration of $A$ is the same in both experiments, while the initial concentration of $B$ is varied. The rates of reaction for these two sets of initial conditions are given by 
# 
# <center>$v(0)_{1} = -\frac{1}{\nu_{A}}\Bigr(\frac{d[A]}{dt}\Bigr|_{\substack{t=0}}\Bigr)_{1} = k[A]_{0}^{m_{A}}([B]_{0})_{1}^{m_{B}}$</center><br/>
# 
# <center>$v(0)_{2} = -\frac{1}{\nu_{A}}\Bigr(\frac{d[A]}{dt}\Bigr|_{\substack{t=0}}\Bigr)_{2} = k[A]_{0}^{m_{A}}([B]_{0})_{2}^{m_{B}}$.</center><br/>
# 
# Because $\nu_{A}$, $k$, $[A]_{0}$, ${m_{A}}$ and $m_{B}$ are the same across the two different measurements, one can use the ratio of the initial rates to solve for $m_{B}$.
# 
# <center>$m_{B} = \frac{\ln{\frac{v(0)_{1}}{v(0)_{2}}}}{\ln{\frac{([B]_{0})_{1}}{([B]_{0})_{2}}}}$</center><br/>
# 
# Another set of experiments can be carried out holding $[B]_{0}$ fixed while varying $[A]_{0}$ to determine $m_{A}$.
# 
# Once the orders are determined, the value of the rate constant $k$ can also be determined through linear regression.
# <center>$\hat{k} = \frac{x^Ty}{\| x\|^2}$</center><br/>
# <center>$x = \begin{bmatrix}([A]_{0}^{m_{A}}[B]_{0}^{m_{B}})_{1} \\ \vdots \\ ([A]_{0}^{m_{A}}[B]_{0}^{m_{B}})_{l}\end{bmatrix}$,$\qquad y = \begin{bmatrix}v(0)_{1} \\ \vdots \\ v(0)_{l}\end{bmatrix}$
# </center><br/>
# 
# ##### Key assumptions and observations
# 1. In the above example the rate law was *a priori* assumed to be of the form
# <center>$r = k[A]^{m_{A}}[B]^{m_{B}}$.</center><br/>
# 2. By using this method, one also assumes that the reactants can be mixed in any desired proportions and the reaction rate can then be measured.
# 3. The method is not appropriate to study *fast* reactions. That is, if the time required to mix the reactants is long compared with the reaction process itself, the rate law cannot be determined using this method.

# ## Computation
# #### Using the method of initial rates and assuming the rate law is of the form: $r = k[S_{1}]^{m_{S_{1}}}[S_{2}]^{m_{S_{2}}}\dots[S_{n}]^{m_{S_{n}}}$.
# 0. Import required modules/libraries

# In[2]:


import math
import numpy


# 1. Compute an approximate $v(0)$ for a particular reaction species $S_{i}$ after having measured concentrations $[S_{i}]$ over a time period from $t=0$ to $t=t'$. Consider the Taylor series expansion of $[S_{i}]$ at 0 and truncated at the quadratic term:
#     <center>
#     <br/>
#     $[S_{i}] \approx [S_{i}]_{0} + \frac{d[S_{i}]}{dt}\Bigr|_{\substack{t=0}}t + \frac{1}{2}\frac{d^2[S_{i}]}{dt^2}\Bigr|_{\substack{t=0}}t^2$
#     </center>
#     <br/>
#     The polynomial relationship between $[S_{i}]$ is theoretically justifiable. It especially holds true for single-reactant reactions and for experiments where all $[S_{i}]$ are equal. If only the measured concentrations for the first 10% of the reaction are used, the error introduced by the approximation is generally small (0.2% for a first-order reactions, 1.2% for a second-order reaction and 3.3% for a third-order reaction). The linear term in the expansion corresponds to the initial rate and can be obtained by solving least-squares problem: $A\hat{x} = b$
#     <center><br/>$\hat{x} = (A^TA)^{-1}A^Tb$</center><br/>
#     <center>$ A = \begin{bmatrix}
#         1 & t_{1}\\
#         1 & t_{2}\\
#         \vdots & \vdots\\
#         1 & t_{10\%}
#         \end{bmatrix}$,
#     $\qquad\hat{x} = \begin{bmatrix}
#     \frac{d[S_{i}]}{dt}\Bigr|_{\substack{t=0}}\\
#     \frac{d^2[S_{i}]}{dt^2}\Bigr|_{\substack{t=0}}
#     \end{bmatrix}$,
#     $\qquad b = \begin{bmatrix}
#         \Bigr(\frac{[S_{i}] - [S_{i}]_{0}}{t}\Bigr)_{1}\\
#         \Bigr(\frac{[S_{i}] - [S_{i}]_{0}}{t}\Bigr)_{2}\\
#         \vdots\\
#         \Bigr(\frac{[S_{i}] - [S_{i}]_{0}}{t}\Bigr)_{10\%}
#     \end{bmatrix}$
# </center><br/>
#     The initial rate can then be computed from the first component of $\hat{x}$.
#     <center>
#     <br/>
#     $v(0) = \frac{1}{\nu_{S_{i}}}\Bigr\lvert\frac{d[S_{i}]}{dt}\Bigr|_{\substack{t=0}}\Bigr\rvert$
#     </center>
#     <br/>

# In[1]:


def compInitialRate(stoichiometric_coefficient, time_points, concentrations, reaction_percentage = 10.0):
    '''Estimates the initial rate of change in concentration of a reaction species

    Parameters
    ----------
    stoich_coefficient : float
        Stoichiometric coefficient of the species whose concentration was measured.
    time_points : array of floats
        Time measurements in seconds taken during the experiment. The number 0 should be at the beginning of array.
    concentrations : array of floats
        Concentration measurements of the reaction species in units of moles per liter taken at the corresponding
        time points. The first element of the array should be initial concentration.
    reaction_percentage : float
        Percentage of reaction specified to determine what initial subset of measurements shall be used.

    Returns
    -------
    float
       initial rate of reaction in units of moles per (liter * second)
    '''
    if reaction_percentage > 10.0:
        print("The reaction percentage specified will probably lead to significant errors in the approximation of the initial rate. Please keep the value of reaction percentage below or at 10%")
        print("\n")
    # Get the initial concentration.
    initial_concentration = concentrations[0]
    
    # Calculate total variation in concentration and multiply by reaction percentage.
    total_conc_variation = numpy.ptp(concentrations)
    restricted_variation = (reaction_percentage/100.0)*total_conc_variation
    
    # The absolute value of the difference between the the concentration at time t and the initial concentration should be
    # less than or equal to the restricted variation.
    concentration_differences = numpy.array(concentrations) - initial_concentration
    concentration_variations = numpy.absolute(concentration_differences)
    upper_index = next((idx
                        for idx, variation in enumerate(concentration_variations)
                        if variation > restricted_variation
                       ), -1)
    
    # If upper index is 1, there are not enough measurements within the percent of the reaction specified and the best
    # approximation available is the slope of the chord obtained from the first concentration measurement.
    if upper_index <= 1:
        print("There are not enough measurements within the percent of the reaction specified and the best approximation available is the slope of the chord obtained from the first concentration measurement.")
        print("\n")
        initial_rate = abs(concentrations[1]-initial_concentration)/(time_points[1]*stoichiometric_coefficient)
        return initial_rate
    
    # Create matrix A with restricted time points
    A = numpy.vstack((numpy.ones(upper_index-1), time_points[1:upper_index])).T
    print(A)
    # Create vector b with restricted ratios: (concentration at time t - initial concentriation)/t
    b = numpy.divide(concentration_differences[1:upper_index], time_points[1:upper_index])
    
    # Solve least-squares problem
    x = numpy.linalg.lstsq(A, b, rcond=None)[0]
    
    # Calculate initial rate
    initial_rate = abs(x[0])/stoichiometric_coefficient
    
    return initial_rate


# 2. Compute the orders $m_{S_{1}}, m_{S_{2}}, \dots, m_{S_{n}}$. Each $m_{S_{i}}$ can be computed from two experiments in which only the initial concentration of $S_{i}$, $[S_{i}]_{0}$ was varied.

# In[7]:


def compOrder(initial_rate_1, initial_rate_2, initial_conc_1, initial_conc_2):
    '''Computes the order a reaction species from two experiments

    Parameters
    ----------
    initial_rate_1 : float
        Initial rate in moles/(liter*second) from experiment 1.
    initial_rate_2 : float
        Initial rate in moles/(liter*second) from experiment 2.
    initial_conc_1 : float
        Initial concentration (in moles per liter) of the species varied in experiment 1.
    initial_conc_2 : float
        Initial concentration (in moles per liter) of the species varied in experiment 2.

    Returns
    -------
    float
       order of reaction species
    '''
    return math.log(initial_rate_1/initial_rate_2)/math.log(initial_conc_1/initial_conc_2)


# 3. Compute the rate constant $k$. This is done here by solving the least squares problem: $\hat{k}x = y$.
# <center><br/>$\hat{k} = \frac{x^Ty}{\| x\|^2}$</center><br/>
# <center>$x = \begin{bmatrix}([S_{1}]_{0}^{m_{S_{1}}}\dots[S_{n}]_{0}^{m_{S_{n}}})_{1} \\ \vdots \\ ([S_{1}]_{0}^{m_{S_{1}}}\dots[S_{n}]_{0}^{m_{S_{n}}})_{l}\end{bmatrix}$, $\qquad y = \begin{bmatrix} v(0)_{1} \\ \vdots \\ v(0)_{l}\end{bmatrix}$
# </center><br/>

# In[14]:


def compRateConstant(initial_concs, initial_rates, orders):
    '''Computes the rate constant from a number experiments using linear regression

    Parameters
    ----------
    initial_concs : 2D-array of floats
        Initial concentrations (in moles per liter) of each species for each experiment structured as
        [experiment_1, experiment_2, ..., experiment_m], where experiment_i is itself an array of 
        initial concentrations for each species structured as [species_1, species_2, ..., species_n].
        Ordering of the species across experiments must be kept the same.
    initial_rates : array of floats
        Initial rates of change in concentration measured in moles/(liter*second) for each experiment.
        The index of the initial rate obtained from experiment_i in the initial_concs argument must be i.
    orders : array of floats
        Order of each reaction species. The index of the order of a species must be the same as the 
        index of the species in the array of initial concentrations.
    
    Example:
    
    orders = [m_S1, m_S2, ..., mSn]
    
    initial_concs = [[conc_S1, conc_S2, ... , conc_Sn], // experiment 1
                     [conc_S1, conc_S2, ... , conc_Sn], // experiment 2
                                        ...,
                     [conc_S1, conc_S2, ... , conc_Sn]] // experiment m
                     
    initial_rates = [experiment_1_rate, experiment_2_rate, ..., experiment_m_rate]

    Returns
    -------
    float
       rate constant
    '''
    # Create vector x with products of powers of initial concentrations
    x = numpy.prod(numpy.power(initial_concs, orders), axis=1).reshape(-1,1)
    
    # Solve least-squares problem
    k = numpy.linalg.lstsq(x, initial_rates, rcond=None)[0][0]
    
    return k


# ##### Additional functions
# 1. Compute $k$, and the orders $m_{S_{1}}, m_{S_{2}}, \dots, m_{S_{n}}$ from a set of at least $n$ independent experiments by solving the least-squares problem: $A\hat{x} = b$.
# <center><br/>$\hat{x} = (A^TA)^{-1}A^Tb$</center><br/>
#     <center>$A = \begin{bmatrix}
#         1 & (\ln{[S_{1}]_{0}})_{1} & (\ln{[S_{2}]_{0}})_{1} & \dots & (\ln{[S_{n}]_{0}})_{1}\\
#         1 & (\ln{[S_{1}]_{0}})_{2} & (\ln{[S_{2}]_{0}})_{2} & \dots & (\ln{[S_{n}]_{0}})_{2}\\
#         \vdots & \vdots & \vdots & \vdots & \vdots\\
#         1 & (\ln{[S_{1}}]_{0})_{l} & (\ln{[S_{2}]_{0}})_{l} & \dots & (\ln{[S_{n}]_{0}})_{l}
#         \end{bmatrix}$,
#     $\qquad\hat{x} = \begin{bmatrix}
#     \ln{k}\\
#     m_{S_{1}}\\
#     m_{S_{2}}\\
#     \vdots\\
#     m_{S_{n}}
#     \end{bmatrix}$, 
#     $\qquad b = \begin{bmatrix}
#        \ln{v(0)_{1}}\\
#         \ln{v(0)_{2}}\\
#         \vdots\\
#         \ln{v(0)_{l}}
#     \end{bmatrix}$,
#     $l\ge n$
# </center><br/>

# In[ ]:


def compLeastSquaresRateLaw(initial_concs, initial_rates):
    '''Computes the rate constant and n orders from a set of at least n experiments using linear regression

    Parameters
    ----------
    initial_concs : 2D-array of floats
        Initial concentrations (in moles per liter) of each species for each experiment structured as
        [experiment_1, experiment_2, ..., experiment_m], where experiment_i is itself an array of 
        initial concentrations for each species structured as [species_1, species_2, ..., species_n].
        Ordering of the species across experiments must be kept the same.
    initial_rates : array of floats
        Initial rates in moles/(liter*second) for each experiment.
        The index of the initial rate obtained from experiment_i in the initial_concs argument must be i.

    Returns
    -------
    array of floats
       [rate constant, m_S1, m_S2, ..., mSn]
    '''
    # Create matrix A with a column of ones and the natural logarithm of each initial concentration.
    A = numpy.vstack((numpy.ones(len(initial_concs)), numpy.log(initial_concs).T)).T

    # Create vector b with the natural logarithm of each initial rate.
    b = numpy.log(initial_rates)

    # Solve the least-squares problem.
    x = numpy.linalg.lstsq(A, b, rcond=None)[0]

    # Replace the first element of x, ln(k), with its antilogarithm to get the value of k.
    x[0] = math.exp(x[0])
    
    return x


# # References
# 1. K. J. Hall, T. I. Quickenden, and D. W. Watts Journal of Chemical Education 1976 53 (8), 493 DOI: 10.1021/ed053p493
# 2. McQuarrie, Donald A. (Donald Allan). Physical Chemistry : a Molecular Approach. Sausalito, Calif. :University Science Books, 1997.
# 
