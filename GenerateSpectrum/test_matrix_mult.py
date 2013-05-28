import numpy as np

nshort = 5
nlong  = 10

x = np.transpose(np.mat(np.ones((nshort))))
y = (np.ones((nshort, nlong)))              # don't need to matricize 2d objects
Nx = len(x)
Ny = len(y)

x[0] = 3
print x
print

print 'test matrix multiplication'
print
print 'x shape, yshape', x.shape, y.shape
print 'x', x
print 'y', y
print

print "Method 1: Simple multiplication and addition; element by element"
z1 = np.zeros(nlong)
for j in range(Nx): z1 += x[j,0]*y[j,:]
print 'z1', z1
print
print

print "Method 2: Abridged version of Method 2; List Comprehension"
z2 = np.zeros(nlong)
z2 = np.sum(np.array( [ x[j,0]*y[j,:] for j in range(nshort) ]), axis=0)
print 'z2', z2
print
print

xnew = x.T
print "Method 3: Matrix Ops", '( shapes' + str(xnew.shape) + ',  ' + str(y.shape)
z3 = np.dot(xnew,y)
print 'z3', z3
print
print
