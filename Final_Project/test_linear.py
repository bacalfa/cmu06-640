#!/usr/bin/env python

'''
Test file for module eosfit
===========================
- Fits data to linear EOS model
- Gets confidence intervals of fitting parameters for linear model
- Plots distribution (histogram) of physical parameters
'''

# Imports
from eosfit import EOS, EOSmodel
import numpy as np

# Data
V = np.array([8., 8.5, 9., 9.6, 10.2, 10.9, 11.6, 12.2, 13., 13.8, 14.5]) # [Ang^3]
E = np.array([-4.65, -5.05, -5.3, -5.48, -5.57, -5.59, -5.575, -5.5, -5.4, -5.3, -5.18]) # [eV/atom]

# Construct EOS model
eos = EOS(V, E, ID='BM4', model=EOSmodel.BM4)
V0_BM4, E0_BM4, B0_BM4 = eos.fit() # No need for initial guess (linear model)
ci_BM4 = eos.get_phys_ci() # Done via simulation (need many samples for good statistical representation)
print '''BM4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
phys ci = {3}
'''.format(V0_BM4, E0_BM4, B0_BM4*160.217, ci_BM4)

# Plot distributions
eos.plot_hist_V0('eosfit_example_BM4_V0_dist.png')
eos.plot_hist_E0('eosfit_example_BM4_E0_dist.png')
eos.plot_hist_B0('eosfit_example_BM4_B0_dist.png')
eos.plot_hist_B0p('eosfit_example_BM4_B0p_dist.png')