import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy.optimize import fsolve

rc('text', usetex=True)


# Define linear system arrays
A = np.array([[1, -3, 9, -27], [1, -1, 1, -1], [1, 1, 1, 1], [1, 2, 4, 8]])
b = np.array([-2, 2, 5, 1])

# Solve system
x = np.linalg.solve(A,b)
print 'Solution of Linear System Ax = b: x = ', x

# Verify solution
bb = np.dot(A,x)
print 'Multiplying A by x gives', bb