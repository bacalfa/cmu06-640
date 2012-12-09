#!/usr/bin/env python

'''
Auxiliary functions for mBM4 EoS model.
Analytical expressions generated by Maple 16
'''

from numpy import sqrt, zeros, asarray, sort

def mBM4_V0(b, c, d, V):
    '''
    Analytical equilibrium volume
    '''
    b = b.astype(complex)
    c = c.astype(complex)
    d = d.astype(complex)
    V0 = zeros((len(b),2))
    
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Ignore complex-to-real casting warning
        V0[:,0] = -(0.3e1*(-c+sqrt((-3*b*d+c**2)))*d-(6*c*d)-0.4e1*(c**2)*(-c+sqrt((-3*b*d+c**2)))/b)/(b**2)
        V0[:,1] = -(-0.3e1*(c+sqrt((-3*b*d+c**2)))*d-(6*c*d)+0.4e1*(c**2)*(c+sqrt((-3*b*d+c**2)))/b)/(b**2)
    V0_vals = zeros((len(b),1))
    
    # Check the volume values that are valid
    for i in range(len(b)):
        v = [V0[i,0], V0[i,1]]
        for vi in v:
            if (vi >= 0.0 and vi >= min(V) and vi <= max(V)):
                V0_vals[i] = vi
                break
    V0_vals_ind = (V0_vals != 0.0).any(axis=1)
    V0_vals = V0_vals[V0_vals_ind]
    return sort(V0_vals_ind), sort(asarray(V0_vals[:,0]))

def mBM4_E0(a, b, c, d, V):
    '''
    Analytical equilibrium energy
    '''
    return sort(a+b*V**(-0.1e1/0.3e1)+c*V**(-0.2e1/0.3e1)+d/V)

def mBM4_B0(b, c, d, V):
    '''
    Analytical equilibrium bulk modulus
    '''
    return sort(0.4e1/0.9e1*b*V**(-0.4e1/0.3e1)+0.10e2/0.9e1*c*V**(-0.5e1/0.3e1)+0.2e1*d/V**2)

def mBM4_B0p(b, c, d, V):
    '''
    Analytical equilibrium first derivative of the bulk modulus
    '''
    return sort(0.2e1/0.3e1*(8*V**(0.4e1/0.3e1)*b+25*V*c+54*d*V**(0.2e1/0.3e1))/V/(V**(0.4e1/0.3e1)*b+4*V*c+9*d*V**(0.2e1/0.3e1)))