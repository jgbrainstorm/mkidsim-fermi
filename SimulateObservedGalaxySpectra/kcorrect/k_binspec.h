#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <kcorrect.h>

/* Create the rmatrix, a lookup table which speeds analysis */
IDL_LONG k_binspec(float *lambda, float *spectrum, float *newlambda, float *newspectrum, IDL_LONG nl, IDL_LONG nnewl);
