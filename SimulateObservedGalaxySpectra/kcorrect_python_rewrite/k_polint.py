#define NRANSI
import numpy
import sys

#void k_polint(float xa[], float ya[], IDL_LONG n, float x, float *y,  float *dy)   # ??? do we care that y and dy are pointers?
def k_polint(xa, ya, n, x, y, dy):
	print "xa", xa, "ya", ya, "n", n, "x", x, "y", y, "dy", dy

	#float *c,*d;                    # ?? do we care that these are pointers?

    ns = long(1)

	dif = numpy.abs(x-xa[1])

	#c=k_vector(1,n);                # ??? is this the same as numpy.arange?
	#d=k_vector(1,n);
	c=numpy.arange(1,n)              # ??? check range is correct
	d=numpy.arange(1,n)
    for i in range(1, n+1):
        if numpy.abs(x- xa[i]) < dif:
			ns  = i
			dif = dift
		#endif
		c[i] = ya[i]
		d[i] = ya[i]

	*y=ya[ns--];                    #??
    for m in range(1,n):
        for i in range(1,(n-m)+1):
			ho = xa[i]  - x
			hp = xa[i+m]- x
			w  = c[i+1] - d[i]
            if (ho- hp) == 0.0 :
				printf("Error in routine polint")
				sys.exit(1)

			den  = w/ den
			d[i] = hp*den
			c[i] = ho*den

		*y += (*dy = ((2* ns) < (n-m) ? c[ns+1] : d[ns--])); #python does not have the ternary operator

	#k_free_vector(d,1,n);
	#k_free_vector(c,1,n);
    del c
    del d

