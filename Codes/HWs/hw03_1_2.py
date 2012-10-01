#!/bin/python

'''
Problem 1.2
===========
'''
from ase.data.molecules import molecule
from jasp import *

# Define molecules
CO = molecule('CO')
CO.set_cell([8, 8, 8], scale_atoms=False)
CO.center()

O2 = molecule('O2')
O2.set_cell([8, 8, 8], scale_atoms=False)
O2.center()

CO2 = molecule('CO2')
CO2.set_cell([8, 8, 8], scale_atoms=False)
CO2.center()

encuts = [250, 350, 450, 500]
DCO, DO2, DCO2, Drxn = [], [], [], []
for encut in encuts:
    with jasp('molecules/hw03/CO-{0}'.format(encut),
        xc='PBE',
        encut=encut,
        ismear=1,
        #sigma=0.01,
        ibrion=2,
        nsw=10,
        atoms=CO) as calc:
        try:
            eCO = CO.get_potential_energy()
            DCO.append(eCO)
        except (VaspSubmitted, VaspQueued):
            eCO = None
    with jasp('molecules/hw03/O2-{0}'.format(encut),
        xc='PBE',
        encut=encut,
        ismear=1,
        #sigma=0.01,
        ibrion=2,
        nsw=10,
        atoms=O2) as calc:
        try:
            eO2 = O2.get_potential_energy()
            DO2.append(eO2)
        except (VaspSubmitted, VaspQueued):
            eO2 = None
    with jasp('molecules/hw03/CO2-{0}'.format(encut),
        xc='PBE',
        encut=encut,
        ismear=1,
        #sigma=0.01,
        ibrion=2,
        nsw=10,
        atoms=CO2) as calc:
        try:
            eCO2 = CO2.get_potential_energy()
            DCO2.append(eCO2)
        except (VaspSubmitted, VaspQueued):
            eCO2 = None
    
    if None not in (eCO, eO2, eCO2):
        dE = eCO2 - eCO - 0.5*eO2
        Drxn.append(dE)

import matplotlib.pyplot as plt
plt.figure()
plt.subplot(221)
plt.plot(encuts, DCO, '-bo')
plt.xlabel('ENCUT (eV)')
plt.ylabel('CO Energy (eV)')
plt.subplot(222)
plt.plot(encuts, DO2, '-ro')
plt.xlabel('ENCUT (eV)')
plt.ylabel('O$_2$ Energy (eV)')
plt.subplot(223)
plt.plot(encuts, DCO2, '-go')
plt.xlabel('ENCUT (eV)')
plt.ylabel('CO$_2$ Energy (eV)')
plt.subplot(224)
plt.plot(encuts, Drxn, '-ko')
plt.xlabel('ENCUT (eV)')
plt.ylabel('CO Oxidation Energy (eV)')
plt.tight_layout()
plt.savefig('images/problem-1-2-convergence.png')