EPS = 1.e-4
JMAX = 14
JMAXP = JMAX + 1
K = 5

#float k_qromo(float (*func)(float),            # !!! Does this get used?
#                       float a, float b,
#							 float (*choose)(float(*)(float), float, float, IDL_LONG))
def k_qromo(lambda, a, b, newlambda_0, newlamba_1 ):

    h = numpy.zeros(JMAXP+1)
    s = numpy.zeros(JMAXP+1)

	h[1]=1.0;
	#for (j=1;j<=JMAX;j++) {
    for j in range(1,JMAX):
        s[j]=k_midpnt(lambda, a,b,j)                                            #
        if j >= K :
			k_polint(h[j-K], s[j-K], K, 0.0, ss, dss)   # output ss, dss        #
            if ( numpy.abs(dss) < (EPS* numpy.abs(ss)) ) or (ss < EPS): return ss
		s[j+1] = s[j]
		h[j+1] = h[j]/ 9.0
	# end for

	print dss,ss
	print "Too many steps in routing qromo"
	return 0.0
}

## fabs = numpy.abs?
