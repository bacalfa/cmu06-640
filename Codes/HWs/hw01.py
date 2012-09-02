'''
06-640: Principles and Applications of Molecular Simulation
Homework 1

Bruno Abreu Calfa
'''

'''
Imports and configurations
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

'''
Problem 4
=========
'''
# Define data
latConst = np.array([3.5 + i*0.05 for i in range(6)])
totE = np.array([-3.649238, -3.696204, -3.719946, -3.723951, -3.711284, -3.68426])

# Fit cubic polynomial
p = np.poly1d(np.polyfit(latConst, totE, 3))

# Get derivative of p and find its roots
dp = np.polyder(p, 1)
z = np.roots(dp)

# Get minimizer of p (second derivative is positive)
dpp = np.polyder(dp, 1)
z_min = z[np.polyval(dpp, z) > 0]
print 'The estimated minimum is = {0:.2f}'.format(z_min[0])

# Plot results
x = np.linspace(latConst[0],latConst[-1])
y = np.polyval(p, x)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(x,y,color='b')
ax.plot(latConst,totE,marker='o',markersize=10,color='r',linestyle='')
ax.plot(z_min,np.polyval(p, z_min),marker='o',markersize=12,color='k')
ax.annotate('Minimum', xy=(z_min, np.polyval(p, z_min) + 0.003), 
            xytext=(z_min - 0.01, np.polyval(p, z_min) + 0.01),
            arrowprops=dict(facecolor='black'),
            )
plt.xlabel('Lattice Constant [$\r{A}$]')
plt.ylabel('Total Energy [eV]')
plt.legend(['Poly Fit', 'Data'], loc='best')
plt.savefig('hw01_prob4.png')
plt.title('Problem 4: Fitting Polynomial to Data')
plt.show()

'''
Problem 5
=========
'''
# Nonlinear function definition
def pb5fcn(x):
    return np.sin(x**2) - 0.5

# Solve the equation
x_0 = 0.5
x_sol = fsolve(pb5fcn, x_0)
print 'One solution of sin(x^2) = 0.5 starting at x_0 = {0} is x* = {1:.4f}'.format(x_0, x_sol[0])

# Plot function and result from fsolve
x = np.linspace(-5,5,1000)
y = pb5fcn(x)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(x,y,color='b')
ax.plot(x_0,pb5fcn(x_0),marker='o',markersize=12,color='r')
ax.annotate('Initial', xy=(x_0, pb5fcn(x_0)), 
            xytext=(x_0 + 1, pb5fcn(x_0)),
            arrowprops=dict(facecolor='red'),
            )
ax.plot(x_sol,pb5fcn(x_sol),marker='o',markersize=12,color='k')
ax.annotate('Root', xy=(x_sol, pb5fcn(x_sol)), 
            xytext=(x_sol + 1, pb5fcn(x_sol)),
            arrowprops=dict(facecolor='black'),
            )
plt.xlabel('x')
plt.ylabel('f(x)')
plt.savefig('hw01_prob5.png')
plt.title('Problem 5: Nonlinear Equation Solution')
plt.show()

'''
Problem 6
=========
'''
# Define linear system arrays
A = np.array([[1, -3, 9, -27], [1, -1, 1, -1], [1, 1, 1, 1], [1, 2, 4, 8]])
b = np.array([-2, 2, 5, 1])

# Solve system
x = np.linalg.solve(A,b)
print 'Solution of Linear System Ax = b: x = ', x

# Verify solution
bb = np.dot(A,x)
print 'Multiplying A by x gives', bb