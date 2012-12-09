#!/usr/bin/env python

'''
Auxiliary functions for mBM5 EoS model.
Analytical expressions generated by Maple 16
'''

from numpy import sqrt, zeros, asarray, sort

def mBM5_V0(b, c, d, e, V):
    '''
    Analytical equilibrium volume
    '''
    b = b.astype(complex)
    c = c.astype(complex)
    d = d.astype(complex)
    e = e.astype(complex)
    V0 = zeros((len(b),3))
    
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Ignore complex-to-real casting warning
        V0[:,0] = -(0.2e1*c*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1-(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1-0.2e1/0.3e1*c/b)**2+0.3e1*d*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1-(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1-0.2e1/0.3e1*c/b)+(4*e))/b
        V0[:,1] = -(2*c*(-(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.6e1)+((9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.6e1)-(0.2e1/0.3e1*c/b)+(0.1e1/0.2e1*1j)*sqrt(0.3e1)*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1+(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1))**2+3*d*(-(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.6e1)+((9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.6e1)-(0.2e1/0.3e1*c/b)+(0.1e1/0.2e1*1j)*sqrt(0.3e1)*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1+(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1))+(4*e))/b
        V0[:,2] = -(2*c*(-(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.6e1)+((9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.6e1)-(0.2e1/0.3e1*c/b)+(-0.1e1/0.2e1*1j)*sqrt(0.3e1)*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1+(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1))**2+3*d*(-(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.6e1)+((9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.6e1)-(0.2e1/0.3e1*c/b)+(-0.1e1/0.2e1*1j)*sqrt(0.3e1)*(0.1e1/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(0.1e1/0.3e1)/0.3e1+(9*d*b-4*c**2)/b*((27*d*c*b)-(54*e*b**2)-(8*c**3)+0.3e1*sqrt(0.3e1)*sqrt((27*d**3*b-9*d**2*c**2-108*d*c*b*e+108*e**2*b**2+32*e*c**3))*b)**(-0.1e1/0.3e1)/0.3e1))+(4*e))/b
    V0_vals = zeros((len(b),1))
    
    # Check the volume values that are valid
    for i in range(len(b)):
        v = [V0[i,0], V0[i,1], V0[i,2]]
        for vi in v:
            if (vi >= 0.0 and vi >= min(V) and vi <= max(V)):
                V0_vals[i] = vi
                break
    V0_vals_ind = (V0_vals != 0.0).any(axis=1)
    V0_vals = V0_vals[V0_vals_ind]
    return sort(V0_vals_ind), sort(asarray(V0_vals[:,0]))

def mBM5_E0(a, b, c, d, e, V):
    '''
    Analytical equilibrium energy
    '''
    return sort(a+b*V**(-0.1e1/0.3e1)+c*V**(-0.2e1/0.3e1)+d/V+e*V**(-0.4e1/0.3e1))

def mBM5_B0(b, c, d, e, V):
    '''
    Analytical equilibrium bulk modulus
    '''
    return sort(0.4e1/0.9e1*b*V**(-0.4e1/0.3e1)+0.10e2/0.9e1*c*V**(-0.5e1/0.3e1)+0.2e1*d/V**2+0.28e2/0.9e1*e*V**(-0.7e1/0.3e1))

def mBM5_B0p(b, c, d, e, V):
    '''
    Analytical equilibrium first derivative of the bulk modulus
    '''
    return sort(0.2e1/0.3e1*(8*b*V+25*c*V**(0.2e1/0.3e1)+54*d*V**(0.1e1/0.3e1)+98*e)/V/(b*V+4*c*V**(0.2e1/0.3e1)+9*d*V**(0.1e1/0.3e1)+16*e))

def mBM5_B0pp(b, c, d, e, V):
    '''
    Analytical equilibrium second derivative of the bulk modulus
    '''
    return sort(4*(2352*e**2+82*c*b*V**(0.5e1/0.3e1)+171*V**(0.4e1/0.3e1)*b*d+294*e*b*V+657*V*c*d+1180*c*V**(0.2e1/0.3e1)*e+2628*V**(0.1e1/0.3e1)*e*d+12*b**2*V**2+150*c**2*V**(0.4e1/0.3e1)+729*d**2*V**(0.2e1/0.3e1))*V**(0.1e1/0.3e1)/(b*V+4*c*V**(0.2e1/0.3e1)+9*d*V**(0.1e1/0.3e1)+16*e)**3)