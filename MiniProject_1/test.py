#!/usr/bin/env python

'''
Test file for module eosfit
'''

# Imports
from eosfit import EOS, EOSmodel
import numpy as np
import matplotlib.pyplot as plt

V = np.array([8., 8.5, 9., 9.6, 10.2, 10.9, 11.6, 12.2, 13., 13.8, 14.5])
E = np.array([-4.65, -5.05, -5.3, -5.48, -5.57, -5.59, -5.575, -5.5, -5.4, -5.3, -5.18])

# MU4 model
eos = EOS(V, E)
p0 = [-6., 2., 5, 10.]
pMU4 = eos.fit(p0) # Initial guess required (nonlinear model)
ciMU4 = eos.get_ci()
EMU4 = eos.MU4(V, pMU4[0], pMU4[1], pMU4[2], pMU4[3])
print pMU4, ciMU4, EMU4

# BM5 model
eos.set_model(EOSmodel.BM5)
pBM5 = eos.fit() # No initial guess required (linear model)
ciBM5 = eos.get_ci()
EBM5 = eos.BM5(V, pBM5[0], pBM5[1], pBM5[2], pBM5[3], pBM5[4])
print pBM5, ciBM5, EBM5

# Plot results
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(V,E,marker='o',markersize=10,color='r',linestyle='')
ax.plot(V,EMU4,color='b')
ax.plot(V,EBM5,color='k')
plt.xlabel('Volume [$\r{A}^3$]')
plt.ylabel('Energy [eV/atom]')
plt.legend(['Data', 'MU4', 'BM5'], loc='best')
plt.savefig('eosfit_example.png')
plt.title('EOSfit Example')
plt.show()