#!/bin/python

'''
Problem 1.1
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

# Create calculators to get the energies
with jasp('molecules/hw03/CO-{0}'.format(250),
    xc='PBE',
    encut=250,
    ismear=1,
    #sigma=0.01,
    ibrion=2,
    nsw=10,
    atoms=CO) as calc:
    try:
        eCO = CO.get_potential_energy()
        print 'CO Forces\n========='
        print CO.get_forces()
    except (VaspSubmitted, VaspQueued):
        eCO = None
with jasp('molecules/hw03/O2-{0}'.format(250),
    xc='PBE',
    encut=250,
    ismear=1,
    #sigma=0.01,
    ibrion=2,
    nsw=10,
    atoms=O2) as calc:
    try:
        eO2 = O2.get_potential_energy()
        print 'O2 Forces\n========='
        print O2.get_forces()
    except (VaspSubmitted, VaspQueued):
        eO2 = None
with jasp('molecules/hw03/CO2-{0}'.format(250),
    xc='PBE',
    encut=250,
    ismear=1,
    #sigma=0.01,
    ibrion=2,
    nsw=10,
    atoms=CO2) as calc:
    try:
        eCO2 = CO2.get_potential_energy()
        print 'CO2 Forces\n========='
        print CO2.get_forces()
    except (VaspSubmitted, VaspQueued):
        eCO2 = None

if None in (eCO, eO2, eCO2):
    pass
else:
    dE = eCO2 - eCO - 0.5*eO2
    print 'Delta E = {0:1.3f} eV'.format(dE)
    print 'Delta E = {0:1.3f} kcal/mol'.format(dE*23.06035)
    print 'Delta E = {0:1.3f} kJ/mol'.format(dE*96.485)
