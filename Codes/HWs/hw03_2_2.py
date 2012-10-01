#!/bin/python

'''
Problem 2.1
===========
'''
import numpy as np
from jasp import *
c = 3e10 # speed of light cm/s
h = 4.135667516e-15 # eV/s

# Get frequencies and calculate ZPE
with jasp('molecules/hw03/CO_vib-{0}'.format(350)) as calc:
    try:
        eCO = calc.get_potential_energy(calc.get_atoms())
        COfreq = calc.get_vibrational_frequencies()
        COZPE = np.sum([0.5*h*f*c for f in COfreq if isinstance(f, float)])
    except (IOError):
        COfreq = None

with jasp('molecules/hw03/O2_vib-{0}'.format(350)) as calc:
    try:
        eO2 = calc.get_potential_energy(calc.get_atoms())
        O2freq = calc.get_vibrational_frequencies()
        O2ZPE = np.sum([0.5*h*f*c for f in O2freq if isinstance(f, float)])
    except (IOError):
        O2freq = None

with jasp('molecules/hw03/CO2_vib-{0}'.format(350)) as calc:
    try:
        eCO2 = calc.get_potential_energy(calc.get_atoms())
        CO2freq = calc.get_vibrational_frequencies()
        CO2ZPE = np.sum([0.5*h*f*c for f in CO2freq if isinstance(f, float)])
    except (IOError):
        CO2freq = None

if None not in (COfreq, O2freq, CO2freq):
    dE = eCO2 - eCO - 0.5*eO2
    dZPE = CO2ZPE - COZPE - 0.5*O2ZPE
    print 'Delta E = {0:1.3f} eV'.format(dE)
    print 'Delta E = {0:1.3f} kcal/mol'.format(dE*23.06035)
    print 'Delta E = {0:1.3f} kJ/mol'.format(dE*96.485)
    print 'Delta ZPE = {0:1.3f} eV'.format(dZPE)
    print 'Delta ZPE = {0:1.3f} kcal/mol'.format(dZPE*23.06035)
    print 'Delta ZPE = {0:1.3f} kJ/mol'.format(dZPE*96.485)
    print 'Delta Total E = {0:1.3f} eV'.format(dE + dZPE)
    print 'Delta Total E = {0:1.3f} kcal/mol'.format((dE + dZPE)*23.06035)
    print 'Delta Total E = {0:1.3f} kJ/mol'.format((dE + dZPE)*96.485)