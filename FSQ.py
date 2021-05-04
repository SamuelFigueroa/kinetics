#!/usr/bin/env python
# coding: utf-8

# # Determining concentrations by absorption spectroscopy
# #### I. Beer's law
# 
# The relation between the concentration of an absorbing analyte and the measured absorbance $A$ or transmittance $T$ of solutions contained in transparent cells that have a path length of $b$ centimeters is described by Beer's law:
# 
# <center>$A = -\log{T} = \log{\frac{P_{0}}{P}} = \epsilon b c$</center>
# <br/>
# 
# where $P_{0}$ is the radiant power in watts incident on the sample and $P$ is the radiant power transmitted by the sample, $\epsilon$ is the molar absorptivity, and $c$ is the concentration of the absorbing substance.
# 
# A mixture containing more than one absorbing substance, i.e. a multicomponent absorbing system, will exhibit a total absorbance,
# 
# <center>
#     $\begin{align} A_{total} &= A_{1} + A_{2} + \dots + A_{n} \\
#     &= \epsilon_{1}bc_1 + \epsilon_{2}bc_2 + \dots + \epsilon_{n}bc_n \end{align}$
# </center>
# 
# #### II. Key assumptions and observations
# 1. Beer's law only applies to mixtures of non-interacting species.
# 2. Significant deviations from the direct (linear) proportionality expressed by the law are often observed under the following conditions:
#     * At high analyte concentrations, usually above $0.01 M$, solute-solvent and solute-solute interactions, or hydrogen bonding can affect the analyte absorptivity. 
#     * At low analyte concentrations but high concentrations of another species, particularly electrolytes.
#     * When changes in the analyte concentration significantly alters the refractive index $n$ of the solution.
#     <center>$A_{observed} = \frac{n}{(n^2+2)^2}\epsilon bc$</center>
#     <br/>
#     * When the analyte associates, dissociates or reacts with the solvent (e.g., an unbuffered solution of a species affected by an acid-base equilibrium) resulting in a multicomponent absorbing system.
#     <center>$A_{observed} = A_{1} + A_{2} + \dots + A_{n}$</center>
#     <br/>
#     * When polychromatic source radiation is employed to make absorbance measurements and the molar absorptivity of the analyte is not constant across the wavelengths considered.
#     <br/>
#     <center>$A_{observed} = 
#     \log{
#     \frac{P_{0\lambda_{1}} + P_{0\lambda_{2}} + \dots + P_{0\lambda_{n}}}
#     {P_{\lambda_{1}} + P_{\lambda_{2}} + \dots + P_{\lambda_{n}}}} 
#     = 
#     \log{
#     \frac{P_{0\lambda_{1}} + P_{0\lambda_{2}} + \dots + P_{0\lambda_{n}}}
#     {P_{0\lambda_{1}}10^{-\epsilon_{\lambda_{1}}bc} + P_{0\lambda_{2}}10^{-\epsilon_{\lambda_{2}}bc} + \dots + P_{0\lambda_{n}}10^{-\epsilon_{\lambda_{n}}bc}}}
#     $</center>
#     <br/>
#     * When stray radiation is present inside the instrument used to measure absorbance.
#     <center>$A_{observed} = \log{\frac{P_{0} + P_{stray}}{P + P_{stray}}}$</center>
#     <br/>
#     * When the cells holding the analyte and blank solutions are not of equal path length and equivalent optical characteristics.
#     <center>$A_{observed} = \epsilon b c + k$</center>
#     <br/>
# 3. For typical organic molecules, strong absorption bands have molar absorptivities ranging from $10^4$ to $10^5M^{-1}cm^{-1}$.
# 
# ## Experimental methods
# #### I. Determining concentrations of absorbing components in a mixture
# In principle, it is possible to determine the concentrations $c_{1}, c_{2}, \dots, c_{n}$ of the individual components of a mixture even if their spectra overlap completely.

# ## Computation
# #### Using principal component analysis on Fourier-transformed spectra
# 0. Import required modules/libraries

# In[9]:


import numpy
import math


# 1. Compute concentrations of standard mixtures. The number of standard mixtures must be equal to or greater than the number of absorbing components. The concentration of each component must be independently varied within the set of standard mixtures to avoid redundancy. To achieve better accuracy, the concentration range of each component across all standard mixtures should closely bracket the range expected in the unknowns. It is also important to maximize the number of components contained in each of the standard mixtures to compensate for the effects of molecular interactions on the measured absorbance.

# 2. Compute the Fourier transform of the absorption spectra divided by the path length (i.e. the optical densities) measured from each of the standard mixtures.

# In[1]:


def compFTOfSpectrum(absorbances, path_length = 1.0):
    '''Computes the Fourier transform of the absorbances measured divided the path length 

    Parameters
    ----------
    absorbances : array of floats
        Array of absorbances in absorbance units from the spectrum measured. 
    path_length : float
        Path length of cuvette holding the sample in centimeters.

    Returns
    -------
    array of floats
       Fourier transform of the optical densities.
    '''
    fourier_transform = numpy.fft.fft(numpy.divide(absorbances, path_length))
    return fourier_transform


# 3. Create a matrix $F$ with $t$ rows and $m$ columns by arranging a subset of consecutive terms $t$ of the computed Fourier transform in a column for each standard mixture $m$. Compute the variance-covariance matrix $A$ ($A = FF^T$) and diagonalize $A$ to obtain the set of $t$ eigenvalues and eigenvectors. Arrange the eigenvectors in order of decreasing eigenvalues.

# In[2]:


def compSortedEigenSet(f):
    '''Computes the set of eigenvalues and eigenvectors of the variance-covariance matrix computed from the input matrix 
    and sorts them in order of decreasing eigenvalue magnitude. 

    Parameters
    ----------
    f : 2D-array of floats
        The input matrix. For example, the Fourier terms of each standard mixture structured as
        // mix_1,   mix_2,  ...,  mix_m 
        [[ f_1_1,   f_1_2,  ...,  f_1_m], // term_1
         [ f_2_1,   f_2_2,  ...,  f_2_m], // term_2
                            ...,
         [ f_t_1,   f_t_2,  ...,  f_t_m], // term_t

    Returns
    -------
    tuple of arrays
       The first element of the tuple is the list of eigenvalues. 
       The second element of the tuple is the 2D numpy array of the corresponding eigenvectors.
    '''
    f = numpy.array(f)
    a = numpy.matmul(f, f.transpose())
    w, v = numpy.linalg.eigh(a)
    eigvals, eigvecs = zip(*sorted(zip(w, v.T), reverse=True))
    return list(eigvals), numpy.array(list(eigvecs)).T


# 4. Compute the value of the indicator function $ind(i)$ at each number $i$ of eigenvectors used. Find the minimum value of this function. The number of eigenvectors to select from the sorted set will be one less than the number of eigenvectors used at this minimum value. Arrange this subset of eigenvectors in a matrix $V$.

# In[ ]:


def compNumEigenVectors(eigenvalues):
    '''Computes the number of eigenvectors to select from the sorted set of eigenvalues by finding 
    the minimum value of the indicator function

    Parameters
    ----------
    eigenvalues : array of floats
        Array of eigenvalues sorted in order of decreasing magnitude. 
    Returns
    -------
    int
       Number of eigenvectors to select from the sorted set.
    '''
    num_eigenvectors = numpy.argmin([eigenvalue/numpy.sum(eigenvalues[i:]) for i, eigenvalue in enumerate(eigenvalues)]) - 1
    return num_eigenvectors


# 5. Compute the matrix $Z$ of projections of each of the vector representations for the standard mixtures onto the orthogonal eigenvectors, $Z = V^TF$.

# 6. Arrange the concentrations of each one of the $n$ components for each of the $m$ standard mixtures in a $n$ by $m$ matrix $C$. Then, compute the proportionality matrix $P = CZ^T(ZZ^T)^{-1}$. Finally, compute the calibration matrix $M = PV^{T}$. The concentrations of each component in an unknown mixture $c_u$ can be obtained by the following relationship after taking the Fourier transform of the mixture's absorbance spectrum $f_u$. 
# <center>$c_u = Mf_u$</center>

# 7. Compute the expected error in the concentrations of the unknown samples, known as the standard error of estimation ($SEE$), by calculating the concentrations of the standard mixtures from their Fourier-transformed spectra using the calibration matrix $M$.
#     <center><br/>$SEE = \frac{\|C_{actual} - C_{estimated} \|_{F}}{\sqrt{n(m-k)}} = \sqrt{\frac{Tr(C^TMF)}{n(m-k)}}$,</center><br/>
#     where $n$ is the number of absorbing components, $m$ the number of standard mixtures and $k$ the number of eigenvectors used.
# 

# In[11]:


def compSEE(concentration_matrix, calibration_matrix, ft_spectra, num_eigenvectors):
    '''Computes the standard error of estimation

    Parameters
    ----------
    concentration_matrix : 2D-array of floats
        Concentrations of each one of the n components for each of the m standard mixtures arranged in a n by m matrix.
    calibration_matrix : 2D-array of floats
        Matrix obtained from the calibration experiments that relates the Fourier transform of a mixture's abosrbance 
        spectrum to the concentrations of each component in an unknown mixture.
    ft_spectra : 2D-array of floats
        Obtained by arranging the subset of t Fourier terms computed from the absorption spectrum
        of each of the m standard mixtures as a t by m matrix.
    num_eigenvectors : int
        Number of eigenvectors used.
    Returns
    -------
    float
       Standard error of estimation (SEE)
    '''
    C = numpy.array(concentration_matrix)
    M = numpy.array(calibration_matrix)
    F = numpy.array(ft_spectra)
    n, m = C.shape
    try:
        SEE = math.sqrt(numpy.trace(numpy.matmul(C.T,numpy.matmul(M,F))).real/(n*(m-num_eigenvectors)))
        return SEE
    except ZeroDivisionError:
        print('''
        The number of eigenvectors is coincidentally the same as the number of standard mixtures.
        This causes a division by zero in the computation.
        SEE cannot be computed in this case.
        ''')
        return math.nan


# # References
# 1. Skoog, D. A., Holler, F. J. and Crouch, S. R., Principles of Instrumental Analysis, 6th Edition, Thomson Brooks/Cole, Belmont, CA, 2007.
# 2. Weismüller, J. A., & Chanady, A. (1992). Quantitative multicomponent analysis of complex mixtures by means of Full Spectrum quantitation and principal component analysis. TrAC Trends in Analytical Chemistry, 11(3), 86-90.
# 3. Brown, C. W., Bump, E. A., & Obremski, R. J. (1986). Accounting for impurities in spectroscopic multicomponent analysis using Fourier vectors. Applied spectroscopy, 40(7), 1023-1031.
# 4. Brown, C. W., Obremski, R. J., & Anderson, P. (1986). Infrared quantitative analysis in the Fourier domain: processing vector representations. Applied spectroscopy, 40(6), 734-742.
# 5. Donahue, S. M., Brown, C. W., & Obremski, R. J. (1988). Multicomponent analysis using Fourier transform infrared and UV spectra. Applied spectroscopy, 42(2), 353-359.
