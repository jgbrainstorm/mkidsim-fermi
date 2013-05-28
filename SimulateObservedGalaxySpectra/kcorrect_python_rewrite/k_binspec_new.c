#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <kcorrect.h>

/*
 * k_binspec.c
 *
 * Bins a spectrum by just integrating over the pixel edges. 
 * Note that this will change a flux density to a flux.
 *
 * Mike Blanton
 * 9/21/2005 */

#define FREEVEC(a) {if((a)!=NULL) free((char *) (a)); (a)=NULL;}

static float *rb_spectrum=NULL;
static float *rb_lambda=NULL;
static IDL_LONG rb_nl;

/* element of projection --> multplication of (redshifted by cr_z) spectrum
 * and the filter */
float rb_pixel(float lambda) 
{
	float sl,spectrum;
	unsigned long i,ip1;

	k_locate(rb_lambda, rb_nl, lambda, &i);
	if(i>=rb_nl-1 || i<0) return(0.);
	ip1=i+1;
	sl=(lambda-rb_lambda[i])/(rb_lambda[ip1]-rb_lambda[i]);

	spectrum=rb_spectrum[i]+sl*(rb_spectrum[ip1]-rb_spectrum[i]);
	
	return(spectrum);
} /* end filter */

/* Create the rmatrix, a lookup table which speeds analysis */
IDL_LONG k_binspec(float *lambda,
									 float *spectrum,
									 float *newlambda,
									 float *newspectrum,
									 IDL_LONG nl,
									 IDL_LONG nnewl)
{
	float lammin,lammax,scale,currlam;
	IDL_LONG i,l,k,v,indxoff;
	
	/* make local copies of vmatrix */
	rb_spectrum=(float *) malloc(nl*sizeof(float));
	rb_lambda=(float *) malloc(nl*sizeof(float));
	rb_nl=nl;
	for(i=0;i<nl;i++) rb_spectrum[i]=spectrum[i];
	for(i=0;i<nl;i++) rb_lambda[i]=lambda[i];
	
	for(i=0;i<nnewl;i++) 
		newspectrum[i]=k_qromo(rb_pixel,newlambda[i], 
													 newlambda[i+1],k_midpnt);
	
	/* deallocate local vmatrix */
	FREEVEC(rb_spectrum);
	FREEVEC(rb_lambda);
	
	return(1);
} /* end create_r */
/* 
 * I have altered this routine to be compatible with zero-offset
 * arrays
 */
void k_locate(float xx[], unsigned long n, float x, unsigned long *j)
{
	unsigned long ju,jm,jl;
	int ascnd;

	jl=-1;
	ju=n;
	ascnd=(xx[n-1] > xx[0]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x > xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	*j=jl;
}
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <kcorrect.h>
#define FUNC(x) ((*func)(x))

float k_midpnt(float (*func)(float), float a, float b, IDL_LONG n)
{
	float x,tnm,sum,del,ddel;
	static float s;
	IDL_LONG it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*(float)tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/(float)tnm)/3.0;
		return s;
	}
}
#undef FUNC
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <kcorrect.h>
#define NRANSI

void k_polint(float xa[], float ya[], IDL_LONG n, float x, float *y, 
							float *dy)
{
	IDL_LONG i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=k_vector(1,n);
	d=k_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				printf("Error in routine polint");
				exit(1);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	k_free_vector(d,1,n);
	k_free_vector(c,1,n);
}
#undef NRANSI
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <kcorrect.h>
#define EPS 1.0e-4
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

float k_qromo(float (*func)(float), float a, float b,
							 float (*choose)(float(*)(float), float, float, IDL_LONG)) 
		newspectrum[i]=k_qromo(rb_pixel,newlambda[i], 
													 newlambda[i+1],k_midpnt);
{
	IDL_LONG j;
	float ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			k_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss) || ss<EPS) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	printf("%le %le\n",dss,ss);
	printf("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
