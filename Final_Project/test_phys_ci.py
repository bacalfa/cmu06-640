#!/usr/bin/env python

'''
Test file for module eosfit
===========================
- Fits data to two EOS models (linear and nonlinear)
- Gets confidence intervals of physical parameters for linear model
- Uses physical parameters estimated from linear model as initial guess for fitting nonlinear model
'''

# Imports
from eosfit import EOS, EOSmodel
import numpy as np

# Data
V = np.array([8., 8.5, 9., 9.6, 10.2, 10.9, 11.6, 12.2, 13., 13.8, 14.5]) # [Ang^3]
E = np.array([-4.65, -5.05, -5.3, -5.48, -5.57, -5.59, -5.575, -5.5, -5.4, -5.3, -5.18]) # [eV/atom]

# Linear EOS model
eos1 = EOS(V, E, ID='mBM4', model=EOSmodel.mBM4)
V0_mBM4, E0_mBM4, B0_mBM4 = eos1.fit() # No need for initial guess (linear model)
ci_mBM4 = eos1.get_phys_ci()
print '''mBM4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
phys ci = {3}
'''.format(V0_mBM4, E0_mBM4, B0_mBM4*160.217, ci_mBM4)

# Nonlinear EOS model
eos2 = EOS(V, E, ID='VI4', model=EOSmodel.VI4)
p0 = [V0_mBM4, E0_mBM4, B0_mBM4, eos1.get_B0p()] # Initial guesses from fitting linear model
V0_VI4, E0_VI4, B0_VI4 = eos2.fit(p0) # Provide custom initial guesses
ci_VI4 = eos2.get_ci()
print '''VI4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
'''.format(V0_VI4, E0_VI4, B0_VI4*160.217, ci_VI4)

# Plot multiple EOS objects
EOS.plotm([eos1, eos2], filename='eosfit_example_mBM4_VI4.png')