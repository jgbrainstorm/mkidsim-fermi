
#float k_midpnt(float (*func)(float), float a, float b, IDL_LONG n)
def k_midpnt(float (*func)(float), a, b, n)

    if n == 1:
		return (b- a)* rb_lambda[0.5* (a+ b)]   # !!! correct function?
    else :
		#for(it=1,j=1;j<n-1;j++) it *= 3;
        it = 1
        for j in range(1,n-1): it *= 3
		tnmi = it
		del  = (b- a)/ (3.0* tnm)
		ddel = del+ del
		x    = a+ 0.5* del
		sum  = 0.0
		#for (j=1;j<=it;j++) {
        for j in range(1,it):
			sum += FUNC(x)
			x   += ddel
			sum += FUNC(x)
			x   += del

		s= (s+ (b- a)* sum/ tnm)/ 3.0
		return s

