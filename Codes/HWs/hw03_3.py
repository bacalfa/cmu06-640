#!/bin/python

'''
Problem 3
=========
'''

from jasp import *
from enthought.mayavi import mlab
from ase.data.vdw import vdw_radii
from ase.data.colors import cpk_colors
from ase.data.molecules import molecule

CO2 = molecule('CO2')
CO2.set_cell([8, 8, 8], scale_atoms=False)
CO2.center()

with jasp('molecules/hw03/CO2-{0}'.format(350)) as calc:
    try:
        atoms = calc.get_atoms()
        x, y, z, cd = calc.get_charge_density()
        mlab.figure(bgcolor=(1, 1, 1))
        
        # plot the atoms as spheres
        for atom in atoms:
            mlab.points3d(atom.x,
                          atom.y,
                          atom.z,
                          scale_factor=vdw_radii[atom.number]/5.,
                          resolution=20,
                          # a tuple is required for the color
                          color=tuple(cpk_colors[atom.number]),
                          scale_mode='none')
        
        # draw the unit cell - there are 8 corners, and 12 connections
        a1, a2, a3 = atoms.get_cell()
        origin = [0, 0, 0]
        cell_matrix = [[origin,  a1],
                       [origin,  a2],
                       [origin,  a3],
                       [a1,      a1 + a2],
                       [a1,      a1 + a3],
                       [a2,      a2 + a1],
                       [a2,      a2 + a3],
                       [a3,      a1 + a3],
                       [a3,      a2 + a3],
                       [a1 + a2, a1 + a2 + a3],
                       [a2 + a3, a1 + a2 + a3],
                       [a1 + a3, a1 + a3 + a2]]
        
        for p1, p2 in cell_matrix:
            mlab.plot3d([p1[0], p2[0]], # x-positions
                        [p1[1], p2[1]], # y-positions
                        [p1[2], p2[2]], # z-positions
                        tube_radius=0.02)
        
        
        # Now plot the charge density
        mlab.contour3d(x, y, z, cd, transparent=True)
        
        # this view was empirically found by iteration
        mlab.view(azimuth=-90, elevation=90, distance='auto')
        
        mlab.savefig('images/co2-centered-cd.png')
        mlab.show()
    except (VaspSubmitted, VaspQueued, IOError):
        pass