## Bins a spectrum by just integrating over the pixel edges.
## Note that this will change a flux density to a flux.
import sys
sys.path.append('kcorrect_python_rewrite/')
import k_binspec_py
import numpy
import pylab
import subprocess
import scipy.integrate

def plotfuncs(wave0,wave1,wave1a, func0, func1,func1a,dwaves=[0,0],
				xmin=None, xmax=None,
				ymin=None, ymax=None,
			  fig_number=0, plot_number=0,  title='kbinspec', base_directory='.', filename_addendum='', show_plot=False):

	fig = pylab.figure(fig_number)
	file_figure_out = base_directory+"/Figures/kbinspec_"+str(plot_number)+"_"+filename_addendum+".png"

	pylab.xlabel('wave')
	pylab.ylabel('func')
	pylab.title(title)
	pylab.xticks(fontsize=10)
	pylab.yticks(fontsize=10)
	if xmin == None : xmin=min(wave0)
	if xmax == None : xmax=max(wave0)
	if ymin == None : ymin=min(func1a)
	if ymax == None : ymax=max(func1a)
	pylab.xlim(xmin=xmin, xmax=xmax)
	pylab.ylim(ymin=ymin, ymax=ymax)

	pylab.plot(wave0,func0,'bo', label="f_sq.; pre-kbinspec;  dw="+str(dwaves[0]))
	pylab.plot(wave1,func1,'r+', label="f_sq.--->f_cube? ; post-kbinspec; dw="+str(dwaves[1]), markersize=20.0)
	pylab.plot(wave1a,func1a,'go', label="f_cube, exact", markersize=10.0)

	pylab.legend()

	pylab.savefig(file_figure_out)
	print "... plot saved to ", file_figure_out

	if show_plot : subprocess.call(["open",file_figure_out])


def func_sq(wave):
	return wave**.5

def func_cube(wave):
	return 0*wave

def test_kbinspec(	nwave = 100, 	dwave0 = 2, 	dwave1 = 5, plot_number=0):

	if dwave1 < dwave0 : bintype = 1
	if dwave1 == dwave0 : bintype = 0
	if dwave1 > dwave0 : bintype = 2
	if dwave1 > dwave0 and dwave1 > 10 : bintype = 3

	

	wave0 = numpy.arange(0, nwave, dwave0); nwave0 = len(wave0)
	wave1 = numpy.arange(0, nwave, dwave1); nwave1 = len(wave1)
	func0 = func_sq(wave0)

	func1 = numpy.zeros(nwave1, dtype='float64')
	k_binspec_py.k_binspec(wave0,func0,wave1,func1)
	

	wave1a = numpy.arange(0, nwave, dwave1); nwave1a = len(wave1a)
	func1a = func_cube(wave1a)

	bintype_name = ['same', 'finer', 'coarser', 'much coarser']
#	plotfuncs(wave0,wave1,wave1a, func0, func1,func1a,
#			plot_number=plot_number, fig_number=plot_number, ymin=-0.1*max(func1a),
#			filename_addendum=bintype_name[bintype],
#			title='kbinspec:'+bintype_name[bintype], dwaves=[dwave0, dwave1], base_directory='.', show_plot=True)
	plotfuncs(wave0,wave1,wave1a, func0,wave0,wave1)

	return


def test_qromb(nwave = 100, 	dwave0 = 1, 	dwave1 = 1, plot_number=0):

	# make wave vector
	wave0 = numpy.arange(0, nwave, dwave0); nwave0 = len(wave0)

	# calc sq func
	func0 = func_sq(wave0)

	# calc cube func
	func0a = func_cube(wave0)

	# integrate to get cube func
	#func0b = numpy.array([scipy.integrate.romberg(func_sq, wave0[ii],wave0[ii+1]) for ii in range(len(wave0)-1)])
	func0b = numpy.array([scipy.integrate.romberg(func_sq, min(wave0),wave0[ii]) for ii in range(nwave-1)])
	#func0b = numpy.array(scipy.integrate.romberg(func_sq, min(wave0),max(wave0)))
	#print len(func0b)

	wave0b = wave0[0:nwave-1]
	print len(wave0), len(wave0b)
	print len(func0), len(func0a), len(func0b)

#	print func0b
	#print min(func0b), max(func0b)
	#print func0a[::-1]
	# compare in plots
	plotfuncs(wave0,wave0,wave0b, func0, func0a,func0b,
			plot_number=plot_number, fig_number=plot_number, ymin=-0.1*max(func0a),
			title='qromb', dwaves=[dwave0, dwave0], base_directory='.', show_plot=True)

	return



def main():
	testing_kbinspec = True
	testing_qromb = False
	if testing_kbinspec:
		test_kbinspec(nwave=100,dwave1=1, plot_number =0)
		test_kbinspec(nwave=100,dwave1=0.1, plot_number =1)
		test_kbinspec(nwave=100,dwave1=10, plot_number =2)
		test_kbinspec(nwave=1000,dwave1=20, plot_number =3)

	if testing_qromb:
		test_qromb(nwave = 100, 	dwave0 = 1, 	dwave1 = 1, plot_number=0)


#-------------------------------
if __name__ == '__main__':
		main()
