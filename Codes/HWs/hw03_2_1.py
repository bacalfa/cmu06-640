#!/bin/python

'''
Problem 2
=========
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

with jasp('molecules/hw03/CO_vib-{0}'.format(350),
    encut=400,
    ismear=0,# Gaussian smearing
    ibrion=6,# finite differences with symmetry
    nfree=2, # central differences (default)
    potim=0.015, # default as well
    ediff=1e-8,
    nsw=1,
    atoms=CO) as calc:
    try:
        CO.get_forces()
        energies, modes = calc.get_vibrational_modes()
        print 'CO Energies\n========'
        for i, e in enumerate(energies):
            print '{0:02d}: {1} eV'.format(i, e)
    except (VaspSubmitted, VaspQueued):
        pass

with jasp('molecules/hw03/O2_vib-{0}'.format(350),
    encut=400,
    ismear=0,# Gaussian smearing
    ibrion=6,# finite differences with symmetry
    nfree=2, # central differences (default)
    potim=0.015, # default as well
    ediff=1e-8,
    nsw=1,
    atoms=O2) as calc:
    try:
        O2.get_forces()
        energies, modes = calc.get_vibrational_modes()
        print 'O2 Energies\n========'
        for i, e in enumerate(energies):
            print '{0:02d}: {1} eV'.format(i, e)
    except (VaspSubmitted, VaspQueued):
        pass
        
with jasp('molecules/hw03/CO2_vib-{0}'.format(350),
    encut=400,
    ismear=0,# Gaussian smearing
    ibrion=6,# finite differences with symmetry
    nfree=2, # central differences (default)
    potim=0.015, # default as well
    ediff=1e-8,
    nsw=1,
    atoms=CO2) as calc:
    try:
        CO2.get_forces()
        energies, modes = calc.get_vibrational_modes()
        print 'CO2 Energies\n========'
        for i, e in enumerate(energies):
            print '{0:02d}: {1} eV'.format(i, e)
    except (VaspSubmitted, VaspQueued):
        pass