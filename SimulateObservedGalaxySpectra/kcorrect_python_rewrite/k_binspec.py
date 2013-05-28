import numpy

#
 #k_binspec.py

 #Bins a spectrum by just integrating over the pixel edges.
 #Note that this will change a flux density to a flux.

 #Mike Blanton
 #9/21/2005 */

static float *rb_spectrum=NULL;
static float *rb_lambda=NULL;
static IDL_LONG rb_nl;

# ============================================================
# element of projection --> multplication of (redshifted by cr_z) spectrum
# and the filter */
# ============================================================
def rb_pixel(lambda):

	k_locate(rb_lambda, rb_nl, lambda, &i)  # returns i, rb_lambda

	if ( i >= (rb_nl- 1)) or (i < 0): return 0.

	ip1 = i+ 1

	sl = (lambda- rb_lambda[i])/ (rb_lambda[ip1]- rb_lambda[i])

	spectrum = rb_spectrum[i]+ sl* (rb_spectrum[ip1]- rb_spectrum[i])

	return spectrum

# END rb_pixel



# ============================================================
# Create the rmatrix, a lookup table which speeds analysis */
# ============================================================

def k_binspec(lambda, spectrum, newlambda, newspectrum, nl, nnewl):
	#float lammin,lammax,scale,currlam;  # these aren't used anywhere ... to be deleted
	#IDL_LONG i,l,k,v,indxoff;  # not used anywhere

	# make local copies of vmatrix
	rb_nl = nl
    rb_spectrum = numpy.array(spectrum)
    rb_lambda   = numpy.array(lambda)

    # integrate
    for i in range(nnewl):
        newspectrum[i] = k_qromo(lambda,newlambda[i], newlambda[i+1],k_midpnt) # lambda should be the argument and the function called

	# deallocate local vmatrix
    del rb_spectrum
    del rb_lambda
	#FREEVEC(rb_spectrum);
	#FREEVEC(rb_lambda);

	return 1
