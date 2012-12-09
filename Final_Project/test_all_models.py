#!/usr/bin/env python

'''
Test file for module eosfit
===========================
- Fits data to all EOS models
- Gets confidence intervals of fitting parameters
- Plots each fit separately and all fits together
'''

# Imports
from eosfit import EOS, EOSmodel
import numpy as np

# Data
V = np.array([8., 8.5, 9., 9.6, 10.2, 10.9, 11.6, 12.2, 13., 13.8, 14.5]) # [Ang^3]
E = np.array([-4.65, -5.05, -5.3, -5.48, -5.57, -5.59, -5.575, -5.5, -5.4, -5.3, -5.18]) # [eV/atom]

plotflag = True

# MU4 model
eos1 = EOS(V, E, ID='MU4')
p0 = [12., -3., 1., 5.] # Order [V0, E0, B0, B0']
V0_MU4, E0_MU4, B0_MU4 = eos1.fit(p0) # Initial guess required (nonlinear model)
ci_MU4 = eos1.get_ci()
E_MU4 = eos1.eval()
R2_MU4 = eos1.get_rsquared()
if plotflag:
    eos1.plot(filename='eosfit_example_MU4.png')
print '''MU4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_MU4, E0_MU4, B0_MU4*160.217, ci_MU4, E_MU4, R2_MU4)

# BM5 model
eos2 = EOS(V, E, ID='BM5', model=EOSmodel.BM5)
V0_BM5, E0_BM5, B0_BM5 = eos2.fit() # No need for initial guess (linear model)
ci_BM5 = eos2.get_ci()
E_BM5 = eos2.eval()
R2_BM5 = eos2.get_rsquared()
if plotflag:
    eos2.plot(filename='eosfit_example_BM5.png')
print '''BM5
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_BM5, E0_BM5, B0_BM5*160.217, ci_BM5, E_BM5, R2_BM5)

# BM4 model
eos3 = EOS(V, E, ID='BM4', model=EOSmodel.BM4)
V0_BM4, E0_BM4, B0_BM4 = eos3.fit() # No need for initial guess (linear model)
ci_BM4 = eos3.get_ci()
E_BM4 = eos3.eval()
R2_BM4 = eos3.get_rsquared()
if plotflag:
    eos3.plot(filename='eosfit_example_BM4.png')
print '''BM4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_BM4, E0_BM4, B0_BM4*160.217, ci_BM4, E_BM4, R2_BM4)

# mBM4 model
eos4 = EOS(V, E, ID='mBM4', model=EOSmodel.mBM4)
V0_mBM4, E0_mBM4, B0_mBM4 = eos4.fit() # No need for initial guess (linear model)
ci_mBM4 = eos4.get_ci()
E_mBM4 = eos4.eval()
R2_mBM4 = eos4.get_rsquared()
if plotflag:
    eos4.plot(filename='eosfit_example_mBM4.png')
print '''mBM4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_mBM4, E0_mBM4, B0_mBM4*160.217, ci_mBM4, E_mBM4, R2_mBM4)

# LOG4 model
eos5 = EOS(V, E, ID='LOG4', model=EOSmodel.LOG4)
V0_LOG4, E0_LOG4, B0_LOG4 = eos5.fit() # No need for initial guess (linear model)
ci_LOG4 = eos5.get_ci()
E_LOG4 = eos5.eval()
R2_LOG4 = eos5.get_rsquared()
if plotflag:
    eos5.plot(filename='eosfit_example_LOG4.png')
print '''LOG4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_LOG4, E0_LOG4, B0_LOG4*160.217, ci_LOG4, E_LOG4, R2_LOG4)

# VI4 model
eos6 = EOS(V, E, ID='VI4', model=EOSmodel.VI4)
p0 = [12., -3., 1., 5.] # Order [V0, E0, B0, B0']
V0_VI4, E0_VI4, B0_VI4 = eos6.fit(p0) # Initial guess required (nonlinear model)
ci_VI4 = eos6.get_ci()
E_VI4 = eos6.eval()
R2_VI4 = eos6.get_rsquared()
if plotflag:
    eos6.plot(filename='eosfit_example_VI4.png')
print '''VI4
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_VI4, E0_VI4, B0_VI4*160.217, ci_VI4, E_VI4, R2_VI4)

# mBM5 model
eos7 = EOS(V, E, ID='mBM5', model=EOSmodel.mBM5)
V0_mBM5, E0_mBM5, B0_mBM5 = eos7.fit() # No need for initial guess (linear model)
ci_mBM5 = eos7.get_ci()
E_mBM5 = eos7.eval()
R2_mBM5 = eos7.get_rsquared()
if plotflag:
    eos7.plot(filename='eosfit_example_mBM5.png')
print '''mBM5
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_mBM5, E0_mBM5, B0_mBM5*160.217, ci_mBM5, E_mBM5, R2_mBM5)

# LOG5 model
eos8 = EOS(V, E, ID='LOG5', model=EOSmodel.LOG5)
V0_LOG5, E0_LOG5, B0_LOG5 = eos8.fit() # No need for initial guess (linear model)
ci_LOG5 = eos8.get_ci()
E_LOG5 = eos8.eval()
R2_LOG5 = eos8.get_rsquared()
if plotflag:
    eos8.plot(filename='eosfit_example_LOG5.png')
print '''LOG5
===
V0 = {0} Ang^3
E0 = {1} eV/atom
B0 = {2} GPa
ci = {3}
E = {4} eV/atom
R^2 = {5}
'''.format(V0_LOG5, E0_LOG5, B0_LOG5*160.217, ci_LOG5, E_LOG5, R2_LOG5)

# Plot multiple EOS objects
if plotflag:
    EOS.plotm([eos1, eos2, eos3, eos4, eos5, eos6, eos7, eos8], filename='eosfit_example_all.png')