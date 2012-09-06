'''
06-640: Principles and Applications of Molecular Simulation
Homework 2

Bruno Abreu Calfa
'''
'''
Imports and configurations
'''
import numpy as np
import matplotlib.pyplot as plt
from ase.data.molecules import molecule
from ase.io import write

linesep = '##########################################################################'

'''
Problem 1.1
===========
'''
atoms = molecule('CH3NO2')
masses = atoms.get_masses()

# MW from database
mw_db = np.sum(masses)

# MW "by hand"
mw_bh = 0.0
for i,atom in enumerate(atoms):
    mw_bh += atom.mass

# Print molecular weights
print 'Problem 1.1'
print '==========='
print 'Molecular weight (database) = {0} u'.format(mw_db)
print 'Molecular weight (by hand) = {0} u'.format(mw_bh)
print linesep

'''
Problem 1.2
===========
'''
# COM from database
com_db = atoms.get_center_of_mass()

# COM "by hand"
com_bh = 0.0
for i,atom in enumerate(atoms):
    com_bh += atom.mass*atom.position
com_bh /= mw_bh

# Print centers of mass
print 'Problem 1.2'
print '==========='
print 'Center of mass (database) = {0}'.format(com_db)
print 'Center of mass (by hand) = {0}'.format(com_bh)
print linesep

'''
Problem 1.3
===========
'''
# MOI from database
moi_db = atoms.get_moments_of_inertia()

# MOI "by hand"
# Inertia tensor matrix
I = np.zeros((3,3))
for i,atom in enumerate(atoms):
    dx = atom.position[0] - com_bh[0]
    dy = atom.position[1] - com_bh[1]
    dz = atom.position[2] - com_bh[2]
    I[0,0] += atom.mass*(dy*dy + dz*dz)
    I[1,1] += atom.mass*(dx*dx + dz*dz)
    I[2,2] += atom.mass*(dx*dx + dy*dy)
    
    I[0,1] -= atom.mass*(dx*dy)
    I[0,2] -= atom.mass*(dx*dz)
    I[1,2] -= atom.mass*(dy*dz)
    
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]

(moi_bh, evects) = np.linalg.eig(I)

# Print centers of mass
print 'Problem 1.3'
print '==========='
print 'Moment of inertia (database) = {0}'.format(moi_db)
print 'Moment of inertia (by hand) = {0}'.format(moi_bh)
print linesep

'''
Problem 1.4
===========
'''
# Print bonds lengths
Cind = 0
print 'Problem 1.4'
print '==========='
for i,atom1 in enumerate(atoms):
    if (atom1.symbol == 'C'):
        Cind = i # Store index for C atom
        k = 1
        for j,atom2 in enumerate(atoms):
            if (atom2.symbol == 'H'):
                print 'Distance C - H{0} = {1} Ang'.format(k,atoms.get_distance(i, j))
                k += 1
        break
print linesep

'''
Problem 1.5
===========
'''
# Get indices for N and O atoms
Nind = 0
Oind = np.zeros(2, np.int)
k = 0
for i,atom1 in enumerate(atoms):
    if (atom1.symbol == 'N'):
        Nind = i
    elif (atom1.symbol == 'O'):
        Oind[k] = i
        k += 1
    if (k == 2):
        break

# Compute vectors of difference of positions
a = atoms.positions[Nind] - atoms.positions[Oind[0]]
b = atoms.positions[Nind] - atoms.positions[Oind[1]]

# Compute angle in degrees
theta_rad = np.arccos(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))
theta_deg = theta_rad*180./np.pi

# Print angle
print 'Problem 1.5'
print '==========='
print 'Angle O - N - O = {0} deg'.format(theta_deg)
print linesep

'''
Problem 1.6
===========
'''
# Generate XYZ file
write('nitromethane.xyz', atoms)

# Print file name
print 'Problem 1.6'
print '==========='
print 'File nitromethane.xyz generated'
print linesep

'''
Problem 1.7
===========
'''
# Display nitromethane
atoms.set_cell([10, 11.5, 12.1])
atoms.center(vacuum=13)
write('nitromethane.png', atoms, show_unit_cell=2,rotation='45x,45y,0z')

# Print file name
print 'Problem 1.7'
print '==========='
print 'File nitromethane.png generated'
print linesep