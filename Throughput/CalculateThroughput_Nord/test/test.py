# testing things for throughput.py

# inverted gaussians
import pylab
import math

gaussian = lambda x,mean,sigma: 3*pylab.exp(-(mean-x)**2/sigma)
sigma 	 = lambda fwhm: 	fwhm/ (pylab.sqrt(8. * math.log(2.) ))

aaa = range(0,10,1)
bbb = range(0,10,1)
aaa.pop(0)
bbb.pop(9)
aaa = pylab.array(aaa)
bbb = pylab.array(bbb)
print aaa
print bbb
ccc = bbb-aaa
print ccc

xmin = -5.
xmax = 5.
dx = 1e-2
mean = 1
fwhm = 2.
sig = sigma(fwhm)
#sigma = 3
x = pylab.arange(xmin,xmax,dx)
y = -gaussian(x,mean,sig)

fig0 = pylab.figure(0)
pylab.plot(x, y, linewidth=2.0,color='r',label='gauss1')
pylab.xlim(xmin=min(x), xmax=max(x))
pylab.ylim(ymin=min(y), ymax=max(y))
pylab.xlabel('x')
pylab.ylabel('y')
pylab.grid(True)
pylab.rcParams['legend.loc'] = 'best'
pylab.legend()
pylab.show()
