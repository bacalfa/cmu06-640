#!/usr/bin/env python

# Imports
from sys import exit
import numpy as np
from scipy.optimize import curve_fit, fmin_ncg
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
    ID = ''                     # ID for the EOS object
    model = EOSmodel.MU4        # Default EOS model
    _modelPhysParamsDict = {}   # Dictinary of models (keys) and if they explictly contain physical parameters (values)
    V = np.array([])            # Volume data
    E = np.array([])            # Energy data
    _p = np.array([])           # Parameters mean values
    _phys_plo = np.array([])    # Physical parameters lower values in confidence interval
    _phys_pup = np.array([])    # Physical parameters upper values in confidence interval
    _res = np.array([])         # Residuals of nonlinear regression
    _ci = np.array([])          # Confidence intervals of _p
    _alpha = 0.0                # User's quantile for confidence interval calculation
    _pcov = np.array([])        # Estimated covariances of optimal _p
    _delta = np.array([])       # Standard deviations for Confidence Interval calculation
    _dof = 0                    # Degrees of freedom
    _V0 = 0.0                   # Equilibrium volume
    _V0_vals = np.array([])     # Values of equilibrium volume from simulatiions
    _E0 = 0.0                   # Equilibrium energy
    _E0_vals = np.array([])     # Values of equilibrium energy from simulations
    _B0 = 0.0                   # Equilibrium bulk modulus
    _B0_vals = np.array([])     # Values of equilibrium bulk modulus from simulations
    _B0p = 0.0                  # Equilibrium first derivative of bulk modulus
    _B0p_vals = np.array([])    # Values of equilibrium first derivative of bulk modulus from simulations
    _B0pp = 0.0                 # Equilibrium second derivative of bulk modulus
    _B0pp_vals = np.array([])   # Values of equilibrium second derivative of bulk modulus from simulations
    
    def __init__(self, V, E, ID='EOS', model=EOSmodel.MU4):
        '''
        Constructor.
        @param V: Volume data (numpy array)
        @param E: Energy data (numpy array)
        @keyword ID: An ID (tag) for the EOS object (default: 'EOS')
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
        self.V = np.array(V)
        self.E = np.array(E)
        self.ID = ID
        self.model = model
        self._modelPhysParamsDict =    {
                                        EOSmodel.BM5: False,
                                        EOSmodel.BM4: False,
                                        EOSmodel.mBM5: False,
                                        EOSmodel.mBM4: False,
                                        EOSmodel.LOG4: False,
                                        EOSmodel.LOG5: False,
                                        EOSmodel.MU4: True,
                                        EOSmodel.VI4: True,
                                        EOSmodel.MO4: False
                                    }
    
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
    
    def fit(self, p0=None, alpha=0.05):
        '''
        Fits the data to a selected EOS model. If model is nonlinear and no initial guesses are
        passed, then the data are fit to BM4 model first and the result is used as the initial 
        guesses for the selected EOS model.
        @keyword p0: Initial guess for parameters (default: None, for linear models only)
        @keyword alpha: Percentile for models without explicit physical parameters (default: 0.05, which means 95% confidence)
        @requires: SciPy package
        '''
        self._alpha = alpha
        # Fit the data to the selected model
        # Linear models: y = Xp, where y is E and X is matrix made of V terms
        # Nonlinear models: use SciPy's curve_fit function
        if self.model == EOSmodel.BM5: # Linear model
            self._fit_BM5()
        
        elif self.model == EOSmodel.mBM5: # Linear model
            self._fit_mBM5()
        
        elif self.model == EOSmodel.BM4: # Linear model
            self._fit_BM4()
        
        elif self.model == EOSmodel.mBM4: # Linear model
            self._fit_mBM4()
        
        elif self.model == EOSmodel.LOG5: # Linear model
            self._fit_LOG5()
            
        elif self.model == EOSmodel.LOG4: # Linear model
            self._fit_LOG4()
        
        elif self.model == EOSmodel.MU4: # Nonlinear model
            if (p0 == None):
                oldModel = self.model
                self.model = EOSmodel.BM4
                self._fit_BM4()
                p0 = np.array([self._V0, self._E0, self._B0, self._B0p])
                self.model = oldModel
            
            self._fit_MU4(p0)    
        
        elif self.model == EOSmodel.VI4: # Nonlinear model
            if (p0 == None):
                oldModel = self.model
                self.model = EOSmodel.BM4
                self._fit_BM4()
                p0 = np.array([self._V0, self._E0, self._B0, self._B0p])
                self.model = oldModel
            
            self._fit_VI4(p0) 
        
        elif self.model == EOSmodel.MO4: # Nonlinear model
            if (p0 == None):
                oldModel = self.model
                self.model = EOSmodel.BM4
                self._fit_BM4()
                a = self._E0 + 9*self._B0*self._V0*(self._B0p - 1.)**(-2)/2.
                b = -9*self._B0*self._V0*(self._B0p - 1.)**(-2)*np.exp(self._B0p - 1.)
                c = 9*self._B0*self._V0*(self._B0p - 1.)**(-2)*np.exp(2*self._B0p - 2.)/2.
                d = (1. - self._B0p)*self._V0**(-1./3.)
                p0 = np.array([a, b, c, d])
                self.model = oldModel
            
            self._fit_MO4(p0) 
        
        return self._V0, self._E0, self._B0
    
    def _fit_BM5(self):
        '''
        Fits the data using BM5 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),5))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = self.V**(-2./3.)
        X[:,2] = self.V**(-4./3.)
        X[:,3] = 1./self.V**2
        X[:,4] = self.V**(-8./3.)
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d, e = self._p
        da, db, dc, dd, de = self._delta
        self._V0, self._E0, self._B0, self._B0p, self._B0pp = self._sim_EOS(a, da, b, db, c, dc, d, dd, e, de)
        
    def _fit_mBM5(self):
        '''
        Fits the data using mBM5 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),5))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = self.V**(-1./3.)
        X[:,2] = self.V**(-2./3.)
        X[:,3] = 1./self.V
        X[:,4] = self.V**(-4./3.)
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d, e = self._p
        da, db, dc, dd, de = self._delta
        self._V0, self._E0, self._B0, self._B0p, self._B0pp = self._sim_EOS(a, da, b, db, c, dc, d, dd, e, de)
    
    def _fit_BM4(self):
        '''
        Fits the data using BM4 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),4))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = self.V**(-2./3.)
        X[:,2] = self.V**(-4./3.)
        X[:,3] = 1./self.V**2
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d = self._p
        da, db, dc, dd = self._delta
        self._V0, self._E0, self._B0, self._B0p = self._sim_EOS(a, da, b, db, c, dc, d, dd)
    
    def _fit_mBM4(self):
        '''
        Fits the data using mBM4 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),4))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = self.V**(-1./3.)
        X[:,2] = self.V**(-2./3.)
        X[:,3] = 1./self.V
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d = self._p
        da, db, dc, dd = self._delta
        self._V0, self._E0, self._B0, self._B0p = self._sim_EOS(a, da, b, db, c, dc, d, dd)
    
    def _fit_LOG5(self):
        '''
        Fits the data using LOG5 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),5))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = np.log(self.V)
        X[:,2] = np.log(self.V)**2
        X[:,3] = np.log(self.V)**3
        X[:,4] = np.log(self.V)**4
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d, e = self._p
        da, db, dc, dd, de = self._delta
        self._V0, self._E0, self._B0, self._B0p, self._B0pp = self._sim_EOS(a, da, b, db, c, dc, d, dd, e, de)
    
    def _fit_LOG4(self):
        '''
        Fits the data using LOG4 model.
        '''
        # Create regressor matrix (X)
        X = np.zeros((len(self.E),4))
        X[:,0] = np.ones(len(self.E))
        X[:,1] = np.log(self.V)
        X[:,2] = np.log(self.V)**2
        X[:,3] = np.log(self.V)**3
        
        # Compute optimal parameters from pseudo-inverse
        X = np.asmatrix(X)
        self._p = np.linalg.solve(X.T*X, np.asarray(np.dot(X.T, self.E))[0])
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Calculate covariance matrix
        R = np.asmatrix(np.linalg.qr(X, mode='r'))
        Rinv = np.linalg.solve(R, np.eye(np.size(R, axis=0), np.size(R, axis=1)))
        diag_info = np.sum(Rinv*Rinv, axis=1)
        # Store standard deviation
        rmse = np.linalg.norm(self._res)/np.sqrt(self._dof)
        self._delta = np.sqrt(diag_info)*rmse
        a, b, c, d = self._p
        da, db, dc, dd = self._delta
        self._V0, self._E0, self._B0, self._B0p = self._sim_EOS(a, da, b, db, c, dc, d, dd)
    
    def _fit_MU4(self,p0):
        '''
        Fits the data using MU4 model.
        @param p0: Initial guess for the parameters
        '''
        self._p, self._pcov = curve_fit(self._MU4, self.V, self.E, p0)
        self._dof = len(self.E) - len(self._p)
        self._delta = np.sqrt(np.diag(self._pcov))
        self._V0 = self._p[0]
        self._E0 = self._p[1]
        self._B0 = self._p[2]
        self._B0p = self._p[3]
        self._res = self.E - self.eval()
    
    def _fit_VI4(self,p0):
        '''
        Fits the data using VI4 model.
        @param p0: Initial guess for the parameters
        '''
        self._p, self._pcov = curve_fit(self._VI4, self.V, self.E, p0)
        self._dof = len(self.E) - len(self._p)
        self._delta = np.sqrt(np.diag(self._pcov))
        self._V0 = self._p[0]
        self._E0 = self._p[1]
        self._B0 = self._p[2]
        self._B0p = self._p[3]
        self._res = self.E - self.eval()
    
    def _fit_MO4(self,p0):
        '''
        Fits the data using MO4 model.
        @param p0: Initial guess for the parameters
        '''
        self._p, self._pcov = curve_fit(self._MO4, self.V, self.E, p0, maxfev=100000)
        self._dof = len(self.E) - len(self._p)
        self._res = self.E - self.eval()
        # Store standard deviation
        self._delta = np.sqrt(np.diag(self._pcov))
        a, b, c, d = self._p
        da, db, dc, dd = self._delta
        self._V0, self._E0, self._B0, self._B0p = self._sim_EOS(a, da, b, db, c, dc, d, dd)
    
    def eval(self):
        '''
        Evaluates the volume (V) data using the fitted model and returns the energy (E) values.
        '''
        if self.model == EOSmodel.BM5:
            return self._BM5(self.V, self._p[0], self._p[1], self._p[2], self._p[3], self._p[4])
        elif self.model == EOSmodel.mBM5:
            return self._mBM5(self.V, self._p[0], self._p[1], self._p[2], self._p[3], self._p[4])
        elif self.model == EOSmodel.BM4:
            return self._BM4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
        elif self.model == EOSmodel.mBM4:
            return self._mBM4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
        elif self.model == EOSmodel.LOG5:
            return self._LOG5(self.V, self._p[0], self._p[1], self._p[2], self._p[3], self._p[4])
        elif self.model == EOSmodel.LOG4:
            return self._LOG4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
        elif self.model == EOSmodel.MU4:
            return self._MU4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
        elif self.model == EOSmodel.VI4:
            return self._VI4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
        elif self.model == EOSmodel.MO4:
            return self._MO4(self.V, self._p[0], self._p[1], self._p[2], self._p[3])
    
    def get_p(self):
        '''
        Returns the fitting parameters.
        '''
        return self._p
    
    def get_phys_p(self):
        '''
        Returns the physical parameters.
        '''
        if self._modelPhysParamsDict[self.model]:
            return self._p
        else:
            if (len(self._p == 4)):
                return np.array([self._V0, self._E0, self._B0, self._B0p])
            else:
                return np.array([self._V0, self._E0, self._B0, self._B0p, self._B0pp])
    
    def get_ci(self, alpha=0.05):
        '''
        Returns the confidence interval of the estimated parameters.
        You should call 'fit' before calling this method.
        @keyword alpha: Percentile (default: 0.05, which means 95% confidence)
        '''
        self._alpha = alpha
        plo = self._p - t.isf(alpha/2., self._dof)*self._delta
        pup = self._p + t.isf(alpha/2., self._dof)*self._delta
        if self._modelPhysParamsDict[self.model]:
            self._phys_plo = plo
            self._phys_pup = pup
        
        self._ci = np.zeros((len(self._p),2))
        self._ci[:,0] = plo
        self._ci[:,1] = pup
        return self._ci
    
    def get_phys_ci(self, alpha=0.05):
        '''
        Returns the confidence interval of the estimated physical parameters.
        You should call 'fit' before calling this method.
        @keyword alpha: Percentile (default: 0.05, which means 95% confidence)
        '''
        self.ci = self.get_ci(alpha)
        if self._modelPhysParamsDict[self.model]:
            return self.ci
        else:
            phys_p = np.array([self._V0, self._E0, self._B0, self._B0p])
            
            alpha_ind = round(alpha/2.*len(self._V0_vals))
            self._phys_plo = np.array([self._V0_vals[alpha_ind], self._E0_vals[alpha_ind], self._B0_vals[alpha_ind], self._B0p_vals[alpha_ind]])
            self._phys_pup = np.array([self._V0_vals[-alpha_ind], self._E0_vals[-alpha_ind], self._B0_vals[-alpha_ind], self._B0p_vals[-alpha_ind]])
            
            if len(self._delta) == 5:
                np.append(phys_p, [self._B0pp])
                np.append(self._phys_plo, [self._B0pp_vals[alpha_ind]])
                np.append(self._phys_pup, [self._B0pp_vals[-alpha_ind]])
            
            self._phys_plo = phys_p - np.abs(self._phys_plo)
            self._phys_pup = phys_p + np.abs(self._phys_pup)
            
            phys_ci = np.zeros((len(phys_p),2))
            phys_ci[:,0] = self._phys_plo
            phys_ci[:,1] = self._phys_pup
            return phys_ci
    
    def get_rsquared(self):
        '''
        Calculates the coefficient of determination (R^2).
        '''
        err = self.E - self.eval()        
        SSE = np.sum(np.dot(err, err))
        ybar = np.mean(self.E)
        SST = np.sum(np.dot(self.E - ybar, self.E - ybar))
        return 1. - SSE/SST
    
    def get_E0(self):
        '''
        Returns the equilibrium energy from fitted data.
        '''
        return self._E0
    
    def get_V0(self):
        '''
        Returns the equilibrium volume from fitted data.
        '''
        return self._V0
    
    def get_B0(self):
        '''
        Returns the equilibrium bulk modulus from fitted data.
        '''
        return self._B0
    
    def get_B0p(self):
        '''
        Returns the equilibrium first derivative of the bulk modulus from fitted data.
        '''
        return self._B0p
    
    def set_ID(self, ID):
        '''
        Sets the new ID of the EOS object.
        '''
        self.ID = ID
    
    def plot(self,filename=None):
        '''
        Plots the results of the latest fit.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        # Check if lower and upper values for CIs are empty
        if (len(self._phys_plo) < 1):
            # Obtain CIs for physical parameters
            if (self._alpha == 0.0):
                self._alpha = 0.05
            self.get_phys_ci(self._alpha)
        
        # If it won't plot on Linux, uncomment following two lines
        #import matplotlib
        #matplotlib.use('TkAgg')
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        fig.set_size_inches(10,8)
        ax = fig.add_subplot(111)
        ax.plot(self.V,self.E,marker='o',markersize=10,color='r',linestyle='')
        ax.plot(self.V,self.eval(),color='b')
        plt.xlabel('Volume [$\r{A}^3$]')
        plt.ylabel('Energy [eV/atom]')
        plt.legend(['Data', self.ID], loc='best')
        plt.title('''$V_0$ = {0:.4f} ({1:.4f}, {2:.4f}) $\rA^3$
        $E_0$ = {3:.4f} ({4:.4f}, {5:.4f}) eV/atom
        $B_0$ = {6:.4f} ({7:.4f}, {8:.4f}) GPa'''.format(self._V0, self._phys_plo[0], self._phys_pup[0],
                                                         self._E0, self._phys_plo[1], self._phys_pup[1],
                                                         self._B0*160.217, self._phys_plo[2]*160.217, self._phys_pup[2]*160.217))
        if (filename != None):
            plt.savefig(filename)
        plt.show()
    
    @staticmethod
    def plotm(eos_list,filename=None):
        '''
        Plots multiple EOS objects.
        @param eos_list: A list of EOS objects
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        # If it won't plot on Linux, uncomment following two lines
        #import matplotlib
        #matplotlib.use('TkAgg')
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        fig.set_size_inches(10,8)
        ax = fig.add_subplot(111)
        ax.plot(eos_list[0].V,eos_list[0].E,marker='o',markersize=10,color='r',linestyle='')
        plt.hold(b=True)
        legend = ['Data']
        for eos in eos_list:
            ax.plot(eos.V,eos.eval())
            legend.append(eos.ID)
        plt.hold(b=False)
        plt.xlabel('Volume [$\r{A}^3$]')
        plt.ylabel('Energy [eV/atom]')
        plt.legend(legend, loc='best')
        #plt.title('''$V_0$ = {0:.4f} $\pm$ {1:.4f} $\rA^3$
        #$E_0$ = {2:.4f} $\pm$ {3:.4f} eV/atom
        #$B_0$ = {4:.4f} $\pm$ {5:.4f} GPa'''.format(self._V0, self._delta_V0, self._E0, self._delta_E0, self._B0, self._delta_B0))
        if filename != None:
            plt.savefig(filename)
        plt.show()
    
    def plot_hist_V0(self,filename=None):
        '''
        Plots histogram of equilibrium volume values for EOS models without physical fitting paramters.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        plt.hist(self._V0_vals, bins=100, normed=True)
        plt.xlabel('Volume [$\r{A}^3$]')
        plt.ylabel('Normalized Counts')
        plt.title('Histogram of $V_0$')
        if filename != None:
            plt.savefig(filename)
        plt.show()
    
    def plot_hist_E0(self,filename=None):
        '''
        Plots histogram of equilibrium energy values for EOS models without physical fitting paramters.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        plt.hist(self._E0_vals, bins=100, normed=True)
        plt.xlabel('Energy [eV]')
        plt.ylabel('Normalized Counts')
        plt.title('Histogram of $E_0$')
        if filename != None:
            plt.savefig(filename)
        plt.show()
    
    def plot_hist_B0(self,filename=None):
        '''
        Plots histogram of equilibrium bulk modulus values for EOS models without physical fitting paramters.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        plt.hist(160.217*self._B0_vals, bins=100, normed=True)
        plt.xlabel('Bulk Modulus [GPa]')
        plt.ylabel('Normalized Counts')
        plt.title('Histogram of $B_0$')
        if filename != None:
            plt.savefig(filename)
        plt.show()
    
    def plot_hist_B0p(self,filename=None):
        '''
        Plots histogram of equilibrium first derivative of bulk modulus values for EOS models without physical fitting paramters.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        plt.hist(self._B0p_vals, bins=100, normed=True)
        plt.xlabel('First Pressure Derivative of Bulk Modulus [-]')
        plt.ylabel('Normalized Counts')
        plt.title("Histogram of $B_0'$")
        if filename != None:
            plt.savefig(filename)
        plt.show()
    
    def plot_hist_B0pp(self,filename=None):
        '''
        Plots histogram of equilibrium second derivative of bulk modulus values for EOS models without physical fitting paramters.
        @keyword filename: File name (including extension) to save the graph (default: None, no file is saved)
        '''
        import matplotlib.pyplot as plt
        fig = plt.figure(facecolor='white')
        plt.hist(self._B0pp_vals/160.217, bins=100, normed=True)
        plt.xlabel('Second Pressure Derivative of Bulk Modulus [GPa$^{-1}$]')
        plt.ylabel('Normalized Counts')
        plt.title("Histogram of $B_0''$")
        if filename != None:
            plt.savefig(filename)
        plt.show()
        
    def _sim_EOS(self, a, da, b, db, c, dc, d, dd, *args):
        '''
        Simulates the EOS models to obtain average values and standard deviations of physical parameters.
        Only 4- and 5-parameter EOS models are currently supported.
        @param a: Average value of fitting parameter a
        @param da: Standard deviation of fitting parameter a
        @param b: Average value of fitting parameter b
        @param db: Standard deviation of fitting parameter b
        @param c: Average value of fitting parameter c
        @param dc: Standard deviation of fitting parameter c
        @param d: Average value of fitting parameter d
        @param dd: Standard deviation of fitting parameter d
        @param args: Additional parameters and their standard deviations (higher-order models)
        '''
        import random
        random.seed()
        # Sample fitting parameters
        N = 10000 # Number of samples
        t_scores = t.isf(self._alpha/2., self._dof)*self._delta
        a_vals = np.asarray([random.gauss(a, da/t_scores[0]) for i in range(N)])
        b_vals = np.asarray([random.gauss(b, db/t_scores[1]) for i in range(N)])
        c_vals = np.asarray([random.gauss(c, dc/t_scores[2]) for i in range(N)])
        d_vals = np.asarray([random.gauss(d, dd/t_scores[3]) for i in range(N)])
        if (len(args) > 0):
            # Check if len(args) == 2 (5-parameter models)
            if (len(args) == 2):
                e_vals = np.asarray([random.gauss(args[0], args[1]/t_scores[4]) for i in range(N)])
            else:
                print "ERROR: Currently cannot simulate EOS models with more than 5 parameters"
                print "Returning zero values from simulation"
                return [0., 0., 0., 0., [0. for i in range(len(args))]]
        
        if (self.model == EOSmodel.BM4):
            import BM4_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = BM4_aux.BM4_V0(b_vals, c_vals, d_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = BM4_aux.BM4_E0(a_vals, b_vals, c_vals, d_vals, self._V0_vals)
            self._B0_vals = BM4_aux.BM4_B0(b_vals, c_vals, d_vals, self._V0_vals)
            self._B0p_vals = BM4_aux.BM4_B0p(b_vals, c_vals, d_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals)
        elif (self.model == EOSmodel.mBM4):
            import mBM4_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = mBM4_aux.mBM4_V0(b_vals, c_vals, d_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = mBM4_aux.mBM4_E0(a_vals, b_vals, c_vals, d_vals, self._V0_vals)
            self._B0_vals = mBM4_aux.mBM4_B0(b_vals, c_vals, d_vals, self._V0_vals)
            self._B0p_vals = mBM4_aux.mBM4_B0p(b_vals, c_vals, d_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals)
        elif (self.model == EOSmodel.LOG4):
            import LOG4_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = LOG4_aux.LOG4_V0(b_vals, c_vals, d_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = LOG4_aux.LOG4_E0(a_vals, b_vals, c_vals, d_vals, self._V0_vals)
            self._B0_vals = LOG4_aux.LOG4_B0(b_vals, c_vals, d_vals, self._V0_vals)
            self._B0p_vals = LOG4_aux.LOG4_B0p(b_vals, c_vals, d_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals)
        elif (self.model == EOSmodel.MO4):
            import MO4_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = MO4_aux.MO4_V0(b_vals, c_vals, d_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = MO4_aux.MO4_E0(a_vals, b_vals, c_vals, d_vals, self._V0_vals)
            self._B0_vals = MO4_aux.MO4_B0(b_vals, c_vals, d_vals, self._V0_vals)
            self._B0p_vals = MO4_aux.MO4_B0p(b_vals, c_vals, d_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals)
        elif (self.model == EOSmodel.BM5):
            import BM5_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = BM5_aux.BM5_V0(b_vals, c_vals, d_vals, e_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            e_vals = e_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = BM5_aux.BM5_E0(a_vals, b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0_vals = BM5_aux.BM5_B0(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0p_vals = BM5_aux.BM5_B0p(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0pp_vals = BM5_aux.BM5_B0pp(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals), np.mean(self._B0pp_vals)
        elif (self.model == EOSmodel.mBM5):
            import mBM5_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = mBM5_aux.mBM5_V0(b_vals, c_vals, d_vals, e_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            e_vals = e_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = mBM5_aux.mBM5_E0(a_vals, b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0_vals = mBM5_aux.mBM5_B0(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0p_vals = mBM5_aux.mBM5_B0p(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0pp_vals = mBM5_aux.mBM5_B0pp(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals), np.mean(self._B0pp_vals)
        elif (self.model == EOSmodel.LOG5):
            import LOG5_aux
            # Obtain distribution of V0 from exact expressions
            V0_vals_ind, self._V0_vals = LOG5_aux.LOG5_V0(b_vals, c_vals, d_vals, e_vals, self.V)
            
            # Some values of V0 are invalid, so use only the valid values
            a_vals = a_vals[V0_vals_ind]
            b_vals = b_vals[V0_vals_ind]
            c_vals = c_vals[V0_vals_ind]
            d_vals = d_vals[V0_vals_ind]
            e_vals = e_vals[V0_vals_ind]
            
            # Obtain distributions of E0, B0, and B0p
            self._E0_vals = LOG5_aux.LOG5_E0(a_vals, b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0_vals = LOG5_aux.LOG5_B0(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0p_vals = LOG5_aux.LOG5_B0p(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            self._B0pp_vals = LOG5_aux.LOG5_B0pp(b_vals, c_vals, d_vals, e_vals, self._V0_vals)
            return np.mean(self._V0_vals), np.mean(self._E0_vals), np.mean(self._B0_vals), np.mean(self._B0p_vals), np.mean(self._B0pp_vals)
    
    def _BM5(self, V, a, b, c, d, e):
        '''
        Implements the 5-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-2./3.) + c*V**(-4./3.) + d/(V**2) + e*V**(-8./3.)
    
    def _mBM5(self, V, a, b, c, d, e):
        '''
        Implements the Modified 5-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-1./3.) + c*V**(-2./3.) + d/V + e*V**(-4./3.)
    
    def _BM4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-2./3.) + c*V**(-4./3.) + d/(V**2)
    
    def _mBM4(self, V, a, b, c, d):
        '''
        Implements the Modified 4-parameter Birch-Murnaghan EOS model
        '''
        return a + b*V**(-1./3.) + c*V**(-2./3.) + d/V
    
    def _LOG5(self, V, a, b, c, d, e):
        '''
        Implements the 5-parameter Logarithmic EOS model
        '''
        return a + b*np.log(V) + c*(np.log(V))**2 + d*(np.log(V))**3 + e*(np.log(V))**4
    
    def _LOG4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Logarithmic EOS model
        '''
        return a + b*np.log(V) + c*(np.log(V))**2 + d*(np.log(V))**3
    
    def _MU4(self, V, V0, E0, B0, B0p):
        '''
        Implements the 4-parameter Murnaghan EOS model
        '''
        a = E0 - B0*V0/(B0p - 1.)
        return a + B0*V/B0p*(1. + ((V0/V)**B0p)/(B0p - 1.))
    
    def _VI4(self, V, V0, E0, B0, B0p):
        '''
        Implements the 4-parameter Vinet EOS model
        '''
        a = E0 + 4.*B0*V0/(B0p - 1.)**2
        return a - (4.*B0*V0/(B0p - 1.)**2)*(1. - 1.5*(B0p - 1.)*(1. - (V/V0)**(1./3.)))*np.exp(1.5*(B0p - 1.)*(1. - (V/V0)**(1./3.)))
    
    def _MO4(self, V, a, b, c, d):
        '''
        Implements the 4-parameter Morse EOS model
        '''
        return a + b*np.exp(d*V**(1./3.)) + c*np.exp(2.*d*V**(1./3.))