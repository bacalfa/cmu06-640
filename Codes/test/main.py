'''
Created on Aug 28, 2012

@author: bacalfa
'''

#if __name__ == '__main__':
s = ''' ajakjka
lssdks'''
print s

a = [1, 2, 3, 4]
print a, len(a)

b = [4, 3, 2, 1]
c = a + b # Concatenation
d = a*2 # Doubles the list a
print c, d

a = range(4)
b = range(4,10)
c = [1, 'test', 4./6]
print a, b, c

a = [2*i for i in range(4)] # List comprehension
print a

a = (1, 2, 3, 4)
b = list(a)
b[0] = 5
print a, b

print '############################################################'

d = {'key1':23,
     'key2':'test',
     5:[2,3],
     (3,4):'tuple value'}
print d['key1'], d.get(2,None)

print '############################################################'

a = 4
b = 4
if a > b:
    print 'a is greater than b'
elif a >= b:
    print 'a is greater or equal than b'
elif a == b:
    print 'a is equal to b'
elif a <= b:
    print 'a less than or equal to b'
else:
    print 'a is less than b'

print '############################################################'

import myfunc as mf
print mf.square(3)

print '############################################################'

import numpy as np
a = np.array([1,2,3,4])
print a*a # element-wise operation
print np.dot(a,a) # linear-algebra dot product
print len(a)

print '############################################################'

from scipy.optimize import fsolve
def f(x):
    y = 2 - x**2
    return y

x0 = 1.4
x = fsolve(f, x0)
print x
print type(x)

print '############################################################'

import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0,2*np.pi)
y = np.sin(x)

plt.plot(x,y)
plt.plot(x,np.cos(x))
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.legend(['sin(x)', 'cos(x)'], loc='best')
plt.savefig('L02-plot1.png')
plt.show()