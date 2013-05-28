import k_binspec_py
import numpy
import random
import pylab

def plot_spectrum(wavelength, power, iplot, plotname):
    file_figure_out = './Figures/kbinspec_test'+plotname+'.png'
    fig0 = pylab.figure(iplot)
    pylab.plot(wavelength, power, linewidth=1.0,color='k')#,label='gal 1')
    pylab.xlim(xmin=min(wavelength), xmax=max(wavelength))
    pylab.ylim(ymin=min(power), ymax=max(power))
    pylab.xlabel('wavelength [nm or angs]')
    pylab.ylabel('"power" or throughput')
    pylab.title(plotname)
    pylab.grid(True)
    pylab.savefig(file_figure_out)




nbin_old = 27000
#nbin_new = 1000
nreduce = 5


random.seed(99999)

LambdaArr       = numpy.arange(nbin_old,  dtype='float32')
newLambdaArr    = LambdaArr[::nreduce] #numpy.arange(nbin_new,  dtype='float32')
spectrumArr     = numpy.array( [random.randint(1,1000) for x in xrange(nbin_old)],  dtype='float32')
newSpectrumArr  = numpy.zeros(newLambdaArr.shape, dtype='float32')

k_binspec_py.k_binspec(LambdaArr, spectrumArr, newLambdaArr, newSpectrumArr)



print 'LambdaArr', min(LambdaArr), max(LambdaArr)
print 'newLambdaArr', min(newLambdaArr), max(newLambdaArr)
print 'spectrumArr', min(spectrumArr), max(spectrumArr)
print 'newSpectrumArr', min(newSpectrumArr), max(newSpectrumArr)
print
print 'new spec '
#for i in range(len(newSpectrumArr)): print i, newSpectrumArr[i]

iplot = 0
plot_spectrum(LambdaArr, spectrumArr, iplot, 'original'); iplot+=1
plot_spectrum(newLambdaArr, newSpectrumArr, iplot, 'new'); iplot+=1
