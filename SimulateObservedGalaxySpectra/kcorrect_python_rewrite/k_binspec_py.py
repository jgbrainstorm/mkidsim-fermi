import numpy
import sys

# General questions/issues
# 1. there are some pointers.  everything that could be passed by reference received a value.
# 2. is "del" useful in freeing arrays/memory
# 3. does del replace k_free_vector?
# 4. is k_vector(1,n) the same as numpy.arange?
# 5. changed rb_pixel to be explicitly called, rather than passing the function variable between modules
# 6. all variables are passed explicitly between functions, rather than having globals.  will this slow down the code?

# ============================================================
# ============================================================
# ============================================================
#k_binspec.py
#
# Bins a spectrum by just integrating over the pixel edges.
# Note that this will change a flux density to a flux.
#
# Mike Blanton #9/21/2005
# pythonized by Nord
# optimized by Gamper
# ============================================================
# ============================================================
# ============================================================


# ============================================================
# element of projection --> multplication of (redshifted by cr_z) spectrum
# and the filter
# ============================================================
def rb_pixel(lambda_one, rb_lambda, rb_spectrum):

 	if lambda_one < rb_lambda[0] or lambda_one > rb_lambda[len(rb_lambda) - 1]:
 		return 0.
 	elif lambda_one < rb_lambda[rb_pixel.index]:
		rb_pixel.index = k_locate(rb_lambda, len(rb_lambda), lambda_one)
		rb_pixel.factor = (rb_spectrum[rb_pixel.index + 1] - rb_spectrum[rb_pixel.index]) / (rb_lambda[rb_pixel.index + 1] - rb_lambda[rb_pixel.index])
	elif lambda_one <= rb_lambda[rb_pixel.index + 1]:
		pass
	else:
		while lambda_one > rb_lambda[rb_pixel.index + 1]:
			rb_pixel.index += 1
		rb_pixel.factor = (rb_spectrum[rb_pixel.index + 1] - rb_spectrum[rb_pixel.index]) / (rb_lambda[rb_pixel.index + 1] - rb_lambda[rb_pixel.index])
	spectrum = rb_spectrum[rb_pixel.index] + (lambda_one - rb_lambda[rb_pixel.index]) * rb_pixel.factor

# #check if valid
# 	rb_nl = len(rb_lambda)
# 	i = k_locate(rb_lambda, len(rb_lambda), lambda_one)
# 	ip1 = i+ 1
# 	sl = (lambda_one - rb_lambda[i]) / (rb_lambda[ip1] - rb_lambda[i])
# 	spectrumCheck = rb_spectrum[i] + sl * (rb_spectrum[ip1] - rb_spectrum[i])
# 	if i != rb_pixel.index or spectrum != spectrumCheck:
# 		print 'mismatch: ', i, rb_pixel.index, spectrum, spectrumCheck

	#rb_pixel.index = 0
	#rb_pixel.factor = 0

	return spectrum
# END rb_pixel



# ============================================================
# Create the rmatrix, a lookup table which speeds analysis
# ============================================================
def k_binspec(lambda_arr, spectrum, newlambda_arr, newspectrum):

	EPS = 1.e-4
	JMAX = 14
	JMAXP = JMAX + 1
	K = 5

	rb_pixel.index = 0
	rb_pixel.factor = 0

	# compute array lengths
	nl	= len(lambda_arr)
	nnewl = len(newlambda_arr)

	# make local copies of vmatrix
	rb_nl = nl
	rb_spectrum = numpy.array(spectrum)
	rb_lambda   = numpy.array(lambda_arr)

	scale_arr = (rb_spectrum[1:] - rb_spectrum[:-1]) / (rb_lambda[1:] - rb_lambda[:-1])

	binsize = numpy.min(rb_lambda[1:] - rb_lambda[:-1])
	bins = numpy.empty(((rb_lambda[len(rb_lambda) - 1] - rb_lambda[0]) / binsize, 2), numpy.int)
	x = rb_lambda[0]
	l = 0
	r = 0
	for i in range(0, bins.shape[0]):
		while l + 1 < len(rb_lambda) and rb_lambda[l + 1] < x:
			l += 1
		#print r
		while r < len(rb_lambda) and rb_lambda[r] < x + binsize:
			r += 1
			#if r >= len(rb_lambda):
			#	print "err: ", r, "x: ", x, "binsize ",  binsize, "rblamda", rb_lambda[r - 1], "l", l,"x+bin", x + binsize
		bins[i][0] = l
		bins[i][1] = r
		x += binsize

	diff_list = [None, None]
	for i in range(2, K + 1):
		tnm = 3 ** (i - 2)
		diff = numpy.zeros(2 * tnm, dtype=numpy.float64)
		for j in range(1, 2 * tnm):
			diff[j] = diff[j - 1] + (1 + (j & 1)) / (3.0 * tnm);
		diff_list.append(diff)

	# integrate
	for i in range(nnewl - 2):


		if newlambda_arr[i] <= rb_lambda[0] or newlambda_arr[nnewl - 1] > rb_lambda[len(rb_lambda) - 1]:
			continue

		a = newlambda_arr[i]
		b = newlambda_arr[i+1]

		ss = 0.0		# ss and ds need to be instantiated beforehand so they can be passed to k_polint and be updated by k_polint #???
		dss = 0.0

		h = numpy.zeros(JMAXP + 1, dtype=numpy.float64)
		s = numpy.zeros(JMAXP + 1, dtype=numpy.float64)

		h[1] = 1.0
		h[2] = 1 / 9.0
		s[2] = s[1] = (b - a) * rb_pixel(0.5 * (a + b), rb_lambda, rb_spectrum)

		for j in range(2, JMAX + 1):

			tnm = numpy.float64(3 ** (j - 2))

			if len(diff_list) == j:
				diff = numpy.zeros(2 * tnm, dtype=numpy.float64)
				for k in range(1, int(2 * tnm)):
					diff[k] = diff[k - 1] + (1 + (k & 1)) / (3.0 * tnm);
				diff_list.append(diff)

			x = a + (b - a) / (2. * 3. * tnm) + diff_list[j] * (b - a)
			bin_ids = numpy.floor((x - rb_lambda[0]) / binsize).astype(numpy.int)
			bin_index = (bins[bin_ids, 0] + bins[bin_ids, 1]) >> 1
			index_arr = bin_index - (x < rb_lambda[bin_index])
			s[j] = (b - a) / tnm / 3.0 * (rb_spectrum[index_arr] + (x - rb_lambda[index_arr]) * scale_arr[index_arr]).sum()

			if j >= K:
				(ss, dss) = k_polint(h[j - K:], s[j - K:], K, 0.0, ss, dss)
				if (numpy.abs(dss) < (EPS * numpy.abs(ss))) or (ss < EPS):
					newspectrum[i] = ss
					break;

				#print "diff list:", len(diff_list) == j, len(diff_list), j
				#if len(diff_list) == j:
					#diff = numpy.zeros(2 * tnm, dtype=numpy.float64)
					#for k in range(1, 2 * tnm):
						#diff[k] = diff[k - 1] + (1 + (k & 1)) / (3.0 * tnm);
					#diff_list.append(diff)

			s[j + 1] = s[j]
			h[j + 1] = h[j] / 9.0

		if j == JMAX + 1:
			print "Too many steps in routing qromo", j, dss,ss

# #check if valid
# 		newspectrum2 = k_qromo(newlambda_arr[i], newlambda_arr[i+1], rb_lambda, rb_spectrum)
# 		if newspectrum[i] != newspectrum2 and numpy.abs(newspectrum[i] - newspectrum2) > 1e-3:
# 			print 'missmatch: ', newspectrum[i], newspectrum2, newspectrum[i] - newspectrum2

# END k_binspec


# ============================================================
# I have altered this routine to be compatible with zero-offset
# arrays
# ============================================================
def k_locate(xx, n, x): #, j):

	jl = -1
	ju = n
	ascnd = (xx[n-1] > xx[0])			 # !!! logical in python; will it work with the comparison 4 lines below?

	# start while
	while (ju-jl) > 1 :
		jm=(ju+jl) >> 1				 # !!! python has the same operator; does work exactly the same?
		if (x > xx[jm]) == ascnd :
			jl=jm
		else :
			ju=jm
	# end while

	j=jl

	#return							  # OLD: j is output (passed by reference)
	return j							# j cannot be passed by reference as output only; in python it needs to also be input.


# END k_locate


def k_midpnt_sum(x, rb_lambda, rb_spectrum, diff, tnm, i, index_arr, scale_arr):
	for j in range(0, 2 * tnm):
		if x + diff[j] > rb_lambda[i + 1]:
			while x + diff[j] > rb_lambda[i + 1]:
				i += 1
		index_arr[j] = i
	return numpy.sum(rb_spectrum[index_arr] + (x + diff - rb_lambda[index_arr]) * scale_arr[index_arr])

# ============================================================
#float k_midpnt(float (*func)(float), float a, float b, IDL_LONG n)
# ============================================================
def k_midpnt(a, b, n, rb_lambda, rb_spectrum):

	s = 0.
	if n == 1:
		return (b - a) * rb_pixel(0.5 * (a + b), rb_lambda, rb_spectrum)   # !!! rb_pixel explicitly called, rather than passing the function variable between modules
	else :
		tnm 	= 3 ** (n - 2)
		del_val = (b - a) / (3.0 * tnm)

		diff = numpy.zeros(2 * tnm, dtype=numpy.float64)
		index_arr = numpy.empty(2 * tnm, numpy.int)
		for j in range(1, 2 * tnm):
			diff[j] = diff[j - 1] + (1 + (j & 1)) * del_val;
		scale_arr = (rb_spectrum[1:] - rb_spectrum[:-1]) / (rb_lambda[1:] - rb_lambda[:-1])

		#cache this!
		i = k_locate(rb_lambda, len(rb_lambda), a)

		x = a + del_val / 2
		if x < rb_lambda[0] or x > rb_lambda[len(rb_lambda) - 1]:
			sum = 0
		else:
			sum = k_midpnt_sum(x, rb_lambda, rb_spectrum, diff, tnm, i, index_arr, scale_arr)
		s = (b- a) * sum / tnm / 3.0


# #check if valid
# 		s2 = 0.
# 		it = 1
# 		for j in range(1,n-1): it *= 3
# 		tnm 	= it
# 		del_val = (b- a)/ (3.0* tnm)
# 		ddel	= del_val+ del_val
# 		x		= a+ 0.5* del_val
# 		sum  = 0.0
# 		for j in range(1,it+1):
# 			sum += rb_pixel(x, rb_lambda, rb_spectrum)
# 			x   += ddel
# 			sum += rb_pixel(x, rb_lambda, rb_spectrum)
# 			x   += del_val
# 		s2 = (s2 + (b - a) * sum / tnm)/ 3.0
# 		if s != s2 and numpy.abs(s - s2) > 1e-5:
# 			print 'missmatch: ', s, s2, s - s2


		return s

# END k_midpnt



# ============================================================
#void k_polint(float xa[], float ya[], IDL_LONG n, float x, float *y,  float *dy)
# ============================================================

def k_polint(xa, ya, n, x, y, dy): # y and dy are outputs, passed by reference

	ns = long(1)
	#print xa
	dif = numpy.abs(x-xa[1])

	c=numpy.empty(n+2)			 #c=k_vector(1,n);				# ??? is this the same as numpy.arange?
	d=numpy.empty(n+2)
	for i in range(1, n+1):
		dift = numpy.abs(x- xa[i])
		if dift < dif:
			ns  = i
			dif = dift
		#endif
		c[i] = ya[i]
		d[i] = ya[i]
	# end for

	y=ya[ns]
	ns -= 1
	for m in range(1,n):
		for i in range(1,n-m+1) :
			ho = xa[i]  - x
			hp = xa[i+m]- x
			w  = c[i+1] - d[i]
			den = ho- hp
			if (den) == 0.0 :
				print("Error in routine polint")
				sys.exit(1)
			# endif

			den  = w/ den
			d[i] = hp*den
			c[i] = ho*den
		# end for

		# *y += (*dy = ((2* ns) < (n-m) ? c[ns+1] : d[ns--]))   # original
		# dy = c[ns+1] if ((2* ns) < (n-m)) else d[ns--]		# need to change ns--
		if  ((2* ns) < (n-m)) :
			dy = c[ns+1]
		else:
			dy = d[ns]
			ns -= 1

		y += dy
	# end for
	return (y, dy)

# END k_polint

# ============================================================
#float k_qromo(float (*func)(float),
#					   float a, float b,
#							 float (*choose)(float(*)(float), float, float, IDL_LONG))
# ============================================================
def k_qromo( a, b, rb_lambda, rb_spectrum ):

	EPS = 1.e-4
	JMAX = 14
	JMAXP = JMAX + 1
	K = 5

	ss = 0.0		# ss and ds need to be instantiated beforehand os they can be passed to k_polint and be updated by k_polint #???
	dss = 0.0

	h = numpy.zeros(JMAXP+1)
	s = numpy.zeros(JMAXP+1)

	i = k_locate(rb_lambda, len(rb_lambda), a)
	scale_arr = (rb_spectrum[1:] - rb_spectrum[:-1]) / (rb_lambda[1:] - rb_lambda[:-1])

	h[1] = 1.0
	h[2] = h[1]/ 9.0
	s[1] = (b - a) * rb_pixel(0.5 * (a + b), rb_lambda, rb_spectrum)
	s[2] = s[1]
	if a >= rb_lambda[0] and b <= rb_lambda[len(rb_lambda) - 1]:
		for j in range(2, JMAX+1):
			tnm 	= 3 ** (j - 2)

			if len(k_qromo.diff_list) == j:
				diff = numpy.zeros(2 * tnm, dtype=numpy.float64)
				for k in range(1, 2 * tnm):
					diff[k] = diff[k - 1] + (1 + (k & 1)) / (3.0 * tnm);
				k_qromo.diff_list.append(diff)
			if len(k_qromo.index_list) == j:
				k_qromo.index_list.append(numpy.empty(2 * tnm, numpy.int))


			del_val = (b - a) / (3.0 * tnm)
			s[j] = (b- a) * k_midpnt_sum(a + del_val / 2, rb_lambda, rb_spectrum, k_qromo.diff_list[j] * (b - a), tnm, i, k_qromo.index_list[j], scale_arr) / tnm / 3.0

# #check if valid
# 			s2 = k_midpnt(a,b,j, rb_lambda, rb_spectrum)
# 			if s[j] != s2 and numpy.abs(s[j] - s2) > 1e-8:
# 				print 'missmatch: ', s[j], s2, s[j] - s2

			if j >= K :
				(ss, dss) = k_polint(h[j-K:], s[j-K:], K, 0.0, ss, dss)							   # ss output (passed by referenced) from k_polint
				if ( numpy.abs(dss) < (EPS* numpy.abs(ss)) ) or (ss < EPS):
					return ss   # final output
			# endif
			s[j+1] = s[j]
			h[j+1] = h[j]/ 9.0
		# end for

	print "Too many steps in routing qromo", j, dss,ss

	return 0.0

	# output is ss, unless the function fails!

k_qromo.diff_list = [None, None]
k_qromo.index_list = [None, None]
# END k_qromo
