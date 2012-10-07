#!/usr/bin/env python

# Imports
from sys import exit
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import t

class EOSmodel:
    '''
    Emulates an 'enum' structure with the EOS models available.
    '''
    BM5 = 0
    BM4 = 1
    mBM5 = 2
    mBM4 = 3
    LOG4 = 4
    LOG5 = 5
    MU4 = 6
    VI4 = 7
    MO4 = 8

class EOS:
    '''
    Defines several Equations of State (EOS) models
    Models available:
    =================
    BM5: 5-parameter Birch-Murnaghan
    BM4: 4-parameter Birch-Murnaghan
    mBM4: Modified 4-parameter Birch-Murnaghan
    mBM5: Modified 5-parameter Birch-Murnaghan
    LOG4: 4-parameter Logarithmic
    LOG5: 5-parameter Logarithmic
    MU4: 4-parameter Murnaghan
    VI4: 4-parameter Vinet
    MO4: 4-parameter Morse
    ------------------------------------------
    Reference:
    Shang, S-L., Wang, Y., Kim, D., Liu, Z-K. (2010) First-principles thermodynamics from phonon and 
    Debye model: Application to Ni and Ni3Al. Computational Materials Science. 47(4): 1040-1048.
    '''
    
    model = EOSmodel.MU4    # Default EOS model
    V = np.array([])        # Volume data
    E = np.array([])        # Energy data
    p = np.array([])        # Parameters mean values
    res = np.array([])      # Residual of nonlinear regression
    ci = np.array([])       # Confidence intervals of p
    pcov = np.array([])     # Estimated covariance of optimal p
    dof = 0                 # Degrees of freedom
    
    def __init__(self, V, E, model=EOSmodel.MU4):
        '''
        Constructor.
        @keyword V: Volume data (numpy array)
        @keyword E: Energy data (numpy array)
        @keyword model: Model name (default: EOSmodel.MU4. Use class EOSmodel for available models.)
        @requires: NumPy package
        '''
        # Validate inputs
        if len(V) == 0 or len(E) == 0:
            print 'ERROR: Vectors V and E must have length greater than zero'
            exit(0)
        if len(V) != len(E):
            print 'ERROR: Vectors V and E must have the same length'
            exit(0)
        self.__validate_model(model)
        
        # Set data and model
        self.V = V
        self.E = E
        self.model = model
    
    def __validate_model(self, model):
        isValidModel = False
        for mod in dir(EOSmodel):
            if (getattr(EOSmodel, mod) == model):
                isValidModel = True
                break
        if not isValidModel:
            print "ERROR: Model '{0}' is not valid".format(model)
            exit(0)
    
    def set_model(self, model):
        '''
        Sets the current EOS model.
        @param model: The current model
        '''
        self.__validate_model(model)
        self.model = model
    
    def fit(self, p0=None):
        '''
        Fits the data to a selected EOS model.
        @keyword p0: Initial guess for parameters (dafult: None, for linear models only)
        @requires: SciPy package
        '''
        # Fit the data to the selected model
        # Linear models: y = Xp, where y is E and X is matrix made of V terms
        # Nonlinear models: use SciPy's curve_fit function
        if self.model == EOSmodel.BM5: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),5))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = self.V**(-2./3.)
            X[:,2] = self.V**(-4./3.)
            X[:,3] = 1./self.V**2
            X[:,4] = self.V**(-8./3.)
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.BM5(self.V, self.p[0], self.p[1], self.p[2], self.p[3], self.p[4])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.mBM5: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),5))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = self.V**(-1./3.)
            X[:,2] = self.V**(-2./3.)
            X[:,3] = 1./self.V
            X[:,4] = self.V**(-4./3.)
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.BM5(self.V, self.p[0], self.p[1], self.p[2], self.p[3], self.p[4])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.BM4: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),4))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = self.V**(-2./3.)
            X[:,2] = self.V**(-4./3.)
            X[:,3] = 1./self.V**2
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.BM4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.mBM4: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),4))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = self.V**(-1./3.)
            X[:,2] = self.V**(-2./3.)
            X[:,3] = 1./self.V
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.mBM4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.LOG5: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),5))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = np.log(self.V)
            X[:,2] = np.log(self.V)**2
            X[:,3] = np.log(self.V)**3
            X[:,4] = np.log(self.V)**4
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.LOG5(self.V, self.p[0], self.p[1], self.p[2], self.p[3], self.p[4])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.LOG4: # Linear model
            # Create regressor matrix (X)
            X = np.zeros((len(self.E),4))
            X[:,0] = np.ones(len(self.E))
            X[:,1] = np.log(self.V)
            X[:,2] = np.log(self.V)**2
            X[:,3] = np.log(self.V)**3
            
            # Compute optimal parameters from pseudo-inverse
            X = np.asmatrix(X)
            self.p = np.asarray(np.dot(np.linalg.inv(X.T*X)*X.T, self.E))[0]
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.LOG4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            s2 = (np.linalg.norm(self.res))**2/self.dof # Estimated residual variance 
            self.pcov = s2*np.linalg.inv(X.T*X)
            return self.p
        
        elif self.model == EOSmodel.MU4: # Nonlinear model
            self.p, self.pcov = curve_fit(self.MU4, self.V, self.E, p0)
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.MU4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            return self.p
        
        elif self.model == EOSmodel.VI4: # Nonlinear model
            self.p, self.pcov = curve_fit(self.VI4, self.V, self.E, p0)
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.VI4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            return self.p
        
        elif self.model == EOSmodel.MO4: # Nonlinear model
            self.p, self.pcov = curve_fit(self.MO4, self.V, self.E, p0)
            self.dof = len(self.E) - len(self.p)
            self.res = self.E - self.MO4(self.V, self.p[0], self.p[1], self.p[2], self.p[3])
            return self.p
    
    def get_ci(self, alpha=0.95):
        '''
        Returns the confidence interval of the estimated parameters.
        You should call 'fit' before calling this method.
        @keyword alpha: Confidence level (default: 95%)
        '''        
        plo = self.p - t.isf(1. - alpha, self.dof)*np.sqrt(np.diag(self.pcov))
        pup = self.p + t.isf(1. - alpha, self.dof)*np.sqrt(np.diag(self.pcov))
        
        self.ci = np.zeros((len(self.p),2))
        self.ci[:,0] = plo
        self.ci[:,1] = pup
        return self.ci
    
    def BM5(self, V, a, b, c, d, e):
        '''
        Implements the 5-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-2./3.) + c*V**(-4./3.) + d/(V**2) + e*V**(-8./3.)
    
    def mBM5(self, V, a, b, c, d, e):
        '''
        Implements the Modified 5-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-1./3.) + c*V**(-2./3.) + d/V + e*V**(-4./3.)
    
    def BM4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-2./3.) + c*V**(-4./3.) + d/(V**2)
    
    def mBM4(self, V, a, b, c, d):
        '''
        Implements the Modified 4-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-1./3.) + c*V**(-2./3.) + d/V
    
    def LOG5(self, V, a, b, c, d, e):
        '''
        Implements the 5-parameter Logarithmic EOS model
        '''
        return a + b*np.log(V) + c*(np.log(V))**2 + d*(np.log(V))**3 + e*(np.log(V))**4
    
    def LOG4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Logarithmic EOS model
        '''
        return a + b*np.log(V) + c*(np.log(V))**2 + d*(np.log(V))**3
    
    def MU4(self, V, E0, B0, B0p, V0):
        '''
        Implements the 4-parameter Murnaghan EOS model
        '''
        a = E0 - B0*V0/(B0p - 1.)
        return a + B0*V/B0p*(1. + ((V0/V)**B0p)/(B0p - 1.))
    
    def VI4(self, V, E0, B0, B0p, V0):
        '''
        Implements the 4-parameter Vinet EOS model
        '''
        a = E0 + 4.*B0*V0/(B0p - 1.)**2
        return a - (4.*B0*V0/(B0p - 1.)**2)*(1. - 1.5*(B0p - 1.)*(1. - (V/V0)**(1./3.)))*np.exp(1.5*(B0p - 1.)*(1. - (V/V0)**(1./3.)))
    
    def MO4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Morse EOS model
        '''
        return a + b*np.exp(d*V**(1./3.)) + c*np.exp(2.*d*V**(1./3.))